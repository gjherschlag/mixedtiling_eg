#!/usr/bin/env python
import math
import string
import ctypes 
import numpy           as np
import scipy.integrate as si
import time
#import scipy.optimize  as numethods
#from vtutool         import *

class universe:
    def __init__(self, filename,vp=None,rnds=123432,Trial=None, pco=1.0, pco_indx=0,tiling=2):
        # set the parameters
        po2           = 1.0
        self.pco      = pco;
        self.pco_indx = pco_indx
        self.p        = parameters(filename,vp,Trial);
        #
        #initialize arrays
        self.pairp  = np.zeros((len(self.p.sitegeoconl),len(self.p.site_prt),len(self.p.site_prt)), dtype=float)
        self.dpairp = np.zeros((len(self.p.sitegeoconl),len(self.p.site_prt),len(self.p.site_prt)), dtype=float)
        # for the system matrix, we want
        #  1) a 1x3x3 array for the bb
        #  2) a 1x3x3 array for the bc
        #  3) a 1x3^6 array for the cccccc
        # So that the total array must be 9+9+3^6
        self.longt= tiling
        n         = self.longt
        m         = len(self.p.site_prt)
        self.mixedtiling  = np.zeros((m*m+m*m+m**n), dtype=float)
        self.dmixedtiling = np.zeros((m*m+m*m+m**n), dtype=float)
        #
        # 0 - all empty
        self.initialize(0)
        # change reaction rates - set manually here to match Temel, Rueter and Schefler parameters.
        # site
        self.p.siterr[np.where(np.all(self.p.sitert==[0,2,0],axis=1))[0][0]] = 7.2*(pco/7.0) #b e  CO
        self.p.siterr[np.where(np.all(self.p.sitert==[1,2,0],axis=1))[0][0]] = 7.2*(pco/7.0) #c e  CO
        self.p.siterr[np.where(np.all(self.p.sitert==[0,0,2],axis=1))[0][0]] = 0.00028       #b CO e
        self.p.siterr[np.where(np.all(self.p.sitert==[1,0,2],axis=1))[0][0]] = 0.092         #c CO e
        # pair
        #print self.p.pairrt
        self.p.pairrr[np.where(np.all(self.p.pairrt==[0,2,2,1,1],axis=1))[0][0]] = 0.97*(po2/1.0) #bb ee OO
        self.p.pairrr[np.where(np.all(self.p.pairrt==[1,2,2,1,1],axis=1))[0][0]] = 0.97*(po2/1.0) #bb ee OO
        self.p.pairrr[np.where(np.all(self.p.pairrt==[2,2,2,1,1],axis=1))[0][0]] = 0.97*(po2/1.0) #bb ee OO
        #
        self.p.pairrr[np.where(np.all(self.p.pairrt==[0,1,1,2,2],axis=1))[0][0]] = 0.0            #bb OO ee
        self.p.pairrr[np.where(np.all(self.p.pairrt==[1,1,1,2,2],axis=1))[0][0]] = 0.00000028     #bb OO ee
        self.p.pairrr[np.where(np.all(self.p.pairrt==[2,1,1,2,2],axis=1))[0][0]] = 0.0            #bb OO ee
        #
        self.p.pairrr[np.where(np.all(self.p.pairrt==[0,0,1,2,2],axis=1))[0][0]] = 0.000000016    #bb COO ee
        self.p.pairrr[np.where(np.all(self.p.pairrt==[0,1,0,2,2],axis=1))[0][0]] = 0.000000016    #bb OCO ee
        self.p.pairrr[np.where(np.all(self.p.pairrt==[1,0,1,2,2],axis=1))[0][0]] = 0.012          #bc COO ee
        self.p.pairrr[np.where(np.all(self.p.pairrt==[1,1,0,2,2],axis=1))[0][0]] = 0.00000520     #bc OCO ee
        self.p.pairrr[np.where(np.all(self.p.pairrt==[2,0,1,2,2],axis=1))[0][0]] = 0.0017         #cc COO ee
        self.p.pairrr[np.where(np.all(self.p.pairrt==[2,1,0,2,2],axis=1))[0][0]] = 0.0017         #cc OCO ee
    def initialize(self,opt=0):
        if   opt==0:
            #start all empty
            self.mixedtiling[0+self.p.site_prt_dict['CO']*3+self.p.site_prt_dict['CO']]  =1.0
            self.mixedtiling[3*3+self.p.site_prt_dict['CO']*3+self.p.site_prt_dict['CO']]=1.0
            ind  = self.p.site_prt_dict['CO']
            indg = 0;
            for i in range(self.longt):
                indg = ind+3*indg
            self.mixedtiling[3*3*2+indg]=1.0
            #print self.mixedtiling[:9]
            #print self.mixedtiling[9:18]
            #print self.mixedtiling[18:]
    def loadProgram(self,fncpp,recompilecpp_huh=True):
        import os
        from subprocess import call
        if(not os.path.exists(fncpp + '.so')) or recompilecpp_huh:
            err = call(["g++","-m64","-shared","-fPIC","-o", fncpp+'.so', fncpp+'.cpp'])
            if err!=0:
                print "compilation of ", fncpp+'.cpp', " has failed in universe.loadProgram(). Exiting..."
                exit(1)
        self.fs_cpp = ctypes.CDLL(fncpp + '.so') 
    def execute_mixed(self):
        Time = self.p.Time;
        dT   = self.p.dT;
        dt   = self.p.dt;
        self.time = np.arange(0,Time+dT*0.5,dT)
        #
        srrl = ctypes.c_int(self.p.sitert.shape[0])
        prrl = ctypes.c_int(self.p.pairrt.shape[0])
        npty = ctypes.c_int(self.pairp.shape[0])
        lngt = ctypes.c_int(self.longt)
        bigd = ctypes.c_int(np.prod(self.mixedtiling.shape))
        nst  = ctypes.c_int(3)
        lenmt= ctypes.c_int(self.longt)
        #
        self.pairp   = self.pairp.flatten()
        self.dpairp  = self.dpairp.flatten()
        pdp  = self.dpairp.ctypes.data_as(ctypes.POINTER(ctypes.c_double));
        ppp  = self.pairp.ctypes.data_as(ctypes.POINTER(ctypes.c_double));
        #
        dmp  = self.dmixedtiling.ctypes.data_as(ctypes.POINTER(ctypes.c_double));
        mtp  = self.mixedtiling.ctypes.data_as(ctypes.POINTER(ctypes.c_double));
        psrr = self.p.siterr.ctypes.data_as(ctypes.POINTER(ctypes.c_double));
        psrt = self.p.sitert.ctypes.data_as(ctypes.POINTER(ctypes.c_int));
        pprr = self.p.pairrr.ctypes.data_as(ctypes.POINTER(ctypes.c_double));
        pprt = self.p.pairrt.ctypes.data_as(ctypes.POINTER(ctypes.c_int));
        ppsc = self.p.pair_scl.ctypes.data_as(ctypes.POINTER(ctypes.c_double));
        pscn = self.p.sitega.ctypes.data_as(ctypes.POINTER(ctypes.c_int));
        pnop = self.p.nopairs.ctypes.data_as(ctypes.POINTER(ctypes.c_int));
        pmtx = self.p.pairmtx.ctypes.data_as(ctypes.POINTER(ctypes.c_double));
        #
        def f_dmt(y,t,univ):
            py = y.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            univ.fs_cpp.constrain_mt(mtp,bigd);
            univ.fs_cpp.constrain_mt(py,bigd);
            univ.fs_cpp.calc_dmt(dmp,py,pdp,ppp,psrt,psrr,pprt,pprr,pnop,pmtx,pscn,lenmt,srrl,prrl,npty,nst)
            return univ.dmixedtiling
        univ=self; ar = (univ,)
        t0  = time.time()
        hmx=0.0; hmn=0.0
        
        ypairs = si.odeint(f_dmt, self.mixedtiling, self.time, args=ar, Dfun=None, h0=dt, hmax=hmx, hmin=hmn, printmessg=1, mxstep=3000)
        t1     = time.time()-t0
        print "time ", t1
        
        newout=True        
        for i in range(ypairs.shape[0]):
            self.output(newout,self.time[i], ypairs[i])
            newout=False#'''
    def output(self,newout,time,outa):
        if newout:
            fp = open(''.join(['results_mixed',str(self.longt),'/pairp_pco_',str(self.pco_indx),'.csv']),'w')
        else:
            fp = open(''.join(['results_mixed',str(self.longt),'/pairp_pco_',str(self.pco_indx),'.csv']),'a')
        for i in range(len(self.mixedtiling)):
            fp.write(str(outa[i]))
            fp.write(" ")
        fp.write(str(time))
        fp.write("\n")
        fp.close()
class parameters:
    def __init__(self, filename,vp,Trial):
        self.set_params(filename,vp,Trial)
    # Setting parameters
    def set_params(self,filename,vp,Trial):
        i = int;
        #Outfile suffix
        self.outsuffix=''
        #Run time parameters
        self.Time       = self.read_param(filename,'Time')
        self.dT         = self.read_param(filename,'dT')
        self.dt         = self.read_param(filename,'dt')
        #Physical Parameters
        self.Temp       = self.read_param(filename,'Temp')
        #Partial pressures
        self.vapor_prt  = self.read_list(filename,'vaporp')
        self.partialpr  = np.zeros(len(self.vapor_prt))
        for i in enumerate(self.vapor_prt):
            self.partialpr[i[0]]=self.read_param(filename,''.join(['p_',i[1][0]]))
        self.vapor_prt_dict = { self.vapor_prt[i][0] : i for i in range(len(self.vapor_prt)) }
        if(vp!=None):
            for key,val in vp.iteritems():
                self.partialpr[self.vapor_prt_dict[key]]=val
                self.outsuffix=''.join([self.outsuffix,"p",key,str(val),'_'])
        #Surface geom
        self.site_ty      =   self.read_list(filename,'sitetypes')
        self.site_ty_dict = { self.site_ty[i][0] : i for i in range(len(self.site_ty)) }
        self.sitegeo      =   self.read_list(filename,'sitegeo')
        self.sitegeoa     =   np.zeros((len(self.sitegeo),len(self.sitegeo[0])), dtype=np.int32)
        for j in range(self.sitegeoa.shape[0]):
            for k in range(self.sitegeoa.shape[1]):
                self.sitegeoa[j,k] = self.site_ty_dict[self.sitegeo[j][k]]
        self.sitegeocon   =   self.buildsitepairs()
        self.sitegeoconl  =   [0]*len(self.sitegeocon)
        for key,value in self.sitegeocon.iteritems():
            temp = key
            self.sitegeoconl[value]=key
        self.sitega       = np.zeros((len(self.sitegeoconl),len(self.sitegeoconl[0])),dtype=np.int32)
        for j in range(self.sitega.shape[0]):
            for k in range(self.sitega.shape[1]):
                self.sitega[j,k] = self.sitegeoconl[j][k]
        self.pair_scl     = self.findpairscls()
        self.nosites      = np.zeros(len(self.site_ty), dtype=np.int32)
        for i in range(self.sitegeoa.shape[0]):
            for j in range(self.sitegeoa.shape[1]):
                self.nosites[self.sitegeoa[i,j]] += 1
        self.nopairs      = np.zeros(len(self.sitegeoconl), dtype=np.int32)
        self.pairmtx      = np.zeros((len(self.sitegeoconl),2*len(self.sitegeoconl)), dtype=float)
        self.findpairsclp()
        #All particle list
        self.all_prt      =   self.read_list(filename,'allp')
        self.all_prt_dict = { self.all_prt[i][0] : i for i in range(len(self.all_prt)) }
        self.prt_info     =   self.read_partinfo(filename,'E:')
        #Site particle list
        self.site_prt     = []
        tmp = self.read_list(filename,'sitep')
        for i in range(len(tmp)):
            self.site_prt.append(tmp[i])
        self.site_prt_dict = { self.site_prt[i][0] : i for i in range(len(self.site_prt)) }
        self.read_reactions(filename,'Rx:')
        if(Trial!=None):
            self.outsuffix=''.join([self.outsuffix,"Trial",string.zfill(str(Trial),3),'_'])
    def read_param(self,filename,param,ptype=float):
        import string
        f = open(filename,"r")
        for line in f:
            input = line.split()
            if(len(input)>1):
                if(input[0]==param):
                    delchars = ''.join(c for c in map(chr, range(256)) if not (c.isalnum() or c==".") )
                    return ptype(input[1].translate(None,delchars))
        print "Error: Value \'%s\' was not found in \'%s\'" % (param,filename)
        exit(1);
    def read_list(self,filename,param,ptype=None):
        import string
        f = open(filename,"r")
        for line in f:
            inpt = line.split()
            if(len(inpt)>1):
                if(inpt[0]==param):
                    ci = len(inpt);
                    fc = False
                    for i in range(1,len(inpt)):
                        if ~fc and inpt[i][0]=='#':
                            fc=True
                            ci = i;
                    if ptype==None:
                        for i in range(1,ci):
                            inpt[i] = inpt[i].split(',')
                        return inpt[1:ci]
        print "Error: Value \'%s\' was not found in \'%s\'" % (param,filename)
        exit(1);
    def buildsitepairs(self):
        outs = {}
        k = 0;
        for i in range(len(self.sitegeo)):
            ip1 = np.mod(i+1,len(self.sitegeo));
            for j in range(len(self.sitegeo[i])):
                jp1 = np.mod(j+1,len(self.sitegeo[i]));
                ijl   = self.site_ty_dict[self.sitegeo[i][j]]
                ipjl  = self.site_ty_dict[self.sitegeo[ip1][j]]
                ijpl  = self.site_ty_dict[self.sitegeo[i][jp1]]
                ipjpl = self.site_ty_dict[self.sitegeo[ip1][jp1]]
                if ijl<ipjl:
                    t1 = [ijl,ipjl]
                else: 
                    t1 = [ipjl,ijl]
                if ijl<ijpl:
                    t2 = [ijl,ijpl]
                else: 
                    t2 = [ijpl,ijl]
                if ijl<ipjpl:
                    t3 = [ijl,ipjpl]
                else: 
                    t3 = [ipjpl,ijl]
                try:
                    outs[t1[0],t1[1]]
                except:
                    outs[t1[0],t1[1]] = k
                    k += 1
                try:
                    outs[t2[0],t2[1]]
                except:
                    outs[t2[0],t2[1]] = k
                    k += 1
                try:
                    outs[t3[0],t3[1]]
                except:
                    outs[t3[0],t3[1]] = k
                    k += 1
        return outs
    def findpairscls(self):
        outs = np.zeros(len(self.sitegeocon),dtype=float)
        for i in range(self.sitegeoa.shape[0]):
            for j in range(self.sitegeoa.shape[1]):
                sy1=self.sitegeoa[i,j]
                for k in range(5,9,2):
                    i1 = (k-(k/3)*3)-1;
                    j1 = (k/3)-1;
                    io = np.mod(i+i1,self.sitegeoa.shape[0])
                    jo = np.mod(j+j1,self.sitegeoa.shape[1])
                    sy2= self.sitegeoa[io,jo]
                    try:
                        ind = self.sitegeocon[(sy1,sy2)]
                        outs[ind]+=1
                    except:
                        pass
                    try:
                        ind = self.sitegeocon[(sy2,sy1)]
                        outs[ind]+=1 
                    except:
                        pass
        outs[:] = outs[:]/2
        return outs
    def findpairsclp(self):
        halflen = len(self.sitegeoconl)
        for i in range(self.sitegeoa.shape[0]):
            for j in range(self.sitegeoa.shape[1]):
                sy1=self.sitegeoa[i,j]
                for i1 in range(-1,2):
                    for j1 in range(-1,2):
                        if(not(i1==0 and j1==0) and (i1==0 or j1==0)):
                            io = np.mod(i+i1,self.sitegeoa.shape[0])
                            jo = np.mod(j+j1,self.sitegeoa.shape[1])
                            sy2= self.sitegeoa[io,jo]
                            if(sy1<=sy2):
                                mns = sy1; mxs = sy2; io1=i;  jo1=j;  io2=io; jo2=jo;
                                ni =  i1; nj = j1;
                            else:
                                mns = sy2; mxs = sy1; io1=io; jo1=jo; io2=i;  jo2=j;
                                ni = -i1; nj = -j1;
                            ind = self.sitegeocon[(mns,mxs)]
                            self.nopairs[ind]+=1
                            #
                            #print i1, j1, mns, mxs
                            #print ind, self.sitegeocon
                            #print self.nopairs
                            #raw_input()
                            #
                            for i2 in range(-1,2):
                                for j2 in range(-1,2):
                                    if(not(i2==0 and j2==0) and (i2==0 or j2==0)):
                                        #site 1 
                                        if(not(i2==ni and j2==nj)):
                                            ioo = np.mod(io1+i2,self.sitegeoa.shape[0])
                                            joo = np.mod(jo1+j2,self.sitegeoa.shape[1])
                                            so  = self.sitegeoa[ioo,joo]
                                            if(mns<=so):
                                                mns1=mns; mxs1=so;
                                            else:
                                                mxs1=mns; mns1=so;
                                            self.pairmtx[ind,self.sitegeocon[(mns1,mxs1)]]+=1
                                        #site 2
                                        if(not(i2==-ni and j2==-nj)):
                                            ioo = np.mod(io2+i2,self.sitegeoa.shape[0])
                                            joo = np.mod(jo2+j2,self.sitegeoa.shape[1])
                                            so  = self.sitegeoa[ioo,joo]
                                            if(mxs<=so):
                                                mns1=mxs; mxs1=so;
                                            else:
                                                mxs1=mxs; mns1=so;
                                            self.pairmtx[ind,halflen+self.sitegeocon[(mns1,mxs1)]]+=1
    def read_partinfo(self,filename,param,ptype=None):
        import string
        n    = len(self.site_ty)+2 #energy levels at sites, atomization energy, molecular weight
        m    = len(self.all_prt)
        prti = np.zeros((m,n))
        prti[:] = None
        #
        nf = True
        f = open(filename,"r")
        for line in f:
            inpt = line.split()
            if(len(inpt)>1):
                if(inpt[0]==param):
                    l = self.all_prt_dict[inpt[1]]
                    for i in range(0,n):
                        if inpt[i+2]!='-':
                            nf = False
                            if inpt[i+2][0]=='-':
                                prti[l,i] = -float(inpt[i+2][1:])
                            else:
                                prti[l,i] = float(inpt[i+2])
        if nf:
            print "Error: Value \'%s\' was not found in \'%s\'; need this for any dynamics to occur" % (param,filename)
            exit(1);
        return prti
    def read_reactions(self,filename,param):
        R           = 8.3144621
        As          = 8
        n           = len(self.site_ty)
        l           = len(self.sitegeocon)
        siterx      = []
        pairrx      = []
        #
        nf = True
        f = open(filename,"r")
        for line in f:
            inpt = line.split()
            if(len(inpt)>1):
                if(inpt[0]==param):
                    i = 0;
                    while i<len(inpt) and not(i=='#'):
                        i+=1;
                    inpt = inpt[:i]
                    fndm = False
                    ist = []
                    fst = []
                    so  = []
                    ae  = []
                    i = 1;
                    # read form inputs
                    while i<len(inpt) and not(fndm):
                        if inpt[i]=='<->':
                            fndm = True
                        elif inpt[i]!='+':
                            ist.append(inpt[i])
                        i+=1
                    if not(fndm):
                        print "Error: Reaction is ill defined - need '<->' to specify reaction exchange"
                        exit(1);
                    fso = False
                    while i<len(inpt) and not(fso):
                        if inpt[i]=='so:':
                            fso = True
                        elif inpt[i]!='+':
                            fst.append(inpt[i])
                        i+=1
                    if not(fso):
                        print "Error: Sticking coeff / preexponential factors are needed"
                        exit(1);
                    fae = False
                    while i<len(inpt) and not(fae):
                        if inpt[i]=='a:':
                            fae = True
                        else:
                            so.append(inpt[i])
                        i+=1
                    if not(fae):
                        print "Error: Activation energy is needed"
                        print "Format expected is Rx: [foo1] + [foo2] <-> [foo3] so: [bar1] [bar2] a: [bar3]"
                        exit(1);
                    if not(len(so)==2):
                        print "Error: Need two so values for forward and backward reaction"
                        exit(1);
                    ae.append(inpt[i])
                    ################
                    '''print "\nChecking re rx prog"
                    print ae
                    print so
                    print ist
                    print fst
                    raw_input()#'''
                    ################
                    #determine if surface-surface (0) interaction, vapor-surface (1), or surface-vapor (2)
                    #also determine if site or pair reaction
                    intrtype = 0
                    isn = 0
                    fsn = 0
                    delE = 0;
                    MWf  = 0
                    MWb  = 0
                    PPf  = 0
                    PPb  = 0
                    ad   = False
                    des  = False
                    for i in range(len(ist)):
                        if(ist[i][0]=='g'):
                            ad   = True
                            ind  = self.all_prt_dict[ist[i][1:]]
                            MWf  = self.prt_info[ind,1+len(self.site_ty)]
                            ind  = self.vapor_prt_dict[ist[i][1:]]
                            PPf  = self.partialpr[ind]
                        elif(ist[i]!='+'):
                            isn+=1
                    if intrtype>1:
                        print "Error: can only adsorb one gas molecule per reaction"
                        exit(1)
                    for i in range(len(fst)):
                        if(fst[i][0]=='g'):
                            des  = True; #desorption
                            ind  = self.all_prt_dict[fst[i][1:]]
                            MWb  = self.prt_info[ind,1+len(self.site_ty)]
                            ind  = self.vapor_prt_dict[fst[i][1:]]
                            PPb  = self.partialpr[ind]
                        elif(fst[i]!='+'):
                            fsn+=1
                    if fsn!=isn:
                        print "Error: sites are unbalanced in reaction"
                        exit(1);
                    #calculate preexponential factors
                    if fsn==2:
                        fracsy = 0.5;
                    elif fsn==1:
                        fracsy = 0.25
                    else:
                        print "Error: can only handle site and pair reactions"
                        exit(1);
                    if ad:
                        c  = float(so[0])*PPf*As*fracsy*133.322368*6.022*1000/np.sqrt(2*np.pi*0.001*MWf*R)
                        pf = lambda T, oc=c: oc/np.sqrt(T)
                    else:
                        c  = float(so[0])*10**(13)
                        pf = lambda T, oc=c: oc
                    if des:
                        c  = float(so[1])*PPb*As*fracsy*133.322368*6.022*1000/np.sqrt(2*np.pi*0.001*MWb*R)
                        pb = lambda T, oc=c: oc/np.sqrt(T)
                    else:
                        c  = float(so[1])*10**(13)
                        pb = lambda T, oc=c: oc
                    cu = 0; #cutoff away from slow reactions
                    if(isn==1):
                        #site reaction
                        for i in range(len(self.site_ty)):
                            delE=0
                            isite = self.site_ty_dict[self.site_ty[i][0]]
                            #find the change in energy due to the reaction
                            #1.) Atomization
                            for il in range(len(ist)):
                                sind = 0;
                                if ist[il][0]=='g':
                                    sind = 1;
                                cmpp = ist[il][sind:]
                                fndo = False
                                for fl in range(len(fst)):
                                    sind = 0;
                                    if fst[fl][0]=='g':
                                        sind = 1;
                                    cmpf = fst[fl][sind:]
                                    if cmpf==cmpp:
                                        fndo = True
                                if fndo==False:
                                    if not(np.isnan(self.prt_info[self.all_prt_dict[cmpp],0])):
                                        delE += self.prt_info[self.all_prt_dict[cmpp],0]
                            for fl in range(len(fst)):
                                sind = 0;
                                if fst[fl][0]=='g':
                                    sind = 1;
                                cmpp = fst[fl][sind:]
                                fndo = False
                                for il in range(len(ist)):
                                    sind = 0;
                                    if ist[il][0]=='g':
                                        sind = 1;
                                    cmpi = ist[il][sind:]
                                    if cmpi==cmpp:
                                        fndo = True
                                if fndo==False:
                                    if not(np.isnan(self.prt_info[self.all_prt_dict[cmpp],0])):
                                        delE -= self.prt_info[self.all_prt_dict[cmpp],0]
                            #print isite, delE
                            #2.) Site energies
                            allowed=True
                            isyst = -1
                            fsyst = -1
                            for il in range(len(ist)):
                                if ist[il][0]!='g':
                                    isyst = self.site_prt_dict[ist[il]]
                                    if not(np.isnan(self.prt_info[self.all_prt_dict[ist[il]],1+i])):
                                        delE -= self.prt_info[self.all_prt_dict[ist[il]],1+i]
                                    else:
                                        allowed=False
                            for fl in range(len(fst)):
                                if fst[fl][0]!='g':
                                    fsyst = self.site_prt_dict[fst[fl]]
                                    if not(np.isnan(self.prt_info[self.all_prt_dict[fst[fl]],1+i])):
                                        delE += self.prt_info[self.all_prt_dict[fst[fl]],1+i]
                                    else:
                                        allowed=False
                            if isyst<0 or fsyst<0:
                                print "Error: initial or final state not found"
                                exit(1)
                            #print "allowed?", allowed
                            if allowed:
                                aeng = float(ae[0])
                                from sympy.functions.special.delta_functions import Heaviside
                                delEf= Heaviside(delE)*abs(delE)+aeng
                                delEb= Heaviside(-delE)*abs(delE)+aeng
                                rf   = lambda T, dE=delEf, pe=pf: pe(T)*np.exp(-float(dE/(T*R*0.001)))
                                rb   = lambda T, dE=delEb, pe=pb: pe(T)*np.exp(-float(dE/(T*R*0.001)))
                                #print isite, delE, aeng,delEf,delEb,rf,rb, pf, pb
                                #print isite, isyst, fsyst                              
                                if(rf(self.Temp)>cu):
                                    siterx.append([isite,isyst,fsyst,rf])
                                if(rb(self.Temp)>cu):
                                    siterx.append([isite,fsyst,isyst,rb])
                    elif(isn==2):
                        #pair reaction
                        for i in range(len(self.sitegeoconl)):
                            for j in range(2):
                                sy  = np.zeros(2)
                                sy[0] = self.sitegeoconl[i][0]
                                sy[1] = self.sitegeoconl[i][1]
                                if j==1:
                                    sy[0] = self.sitegeoconl[i][1]
                                    sy[1] = self.sitegeoconl[i][0]
                                    if sy[0]==sy[1]:
                                        break;
                                for k in range(2):
                                    delE=0
                                    fsy  = np.zeros(2)
                                    fsy[0] = self.sitegeoconl[i][0]
                                    fsy[1] = self.sitegeoconl[i][1]
                                    if k==1:
                                        fsy[0] = self.sitegeoconl[i][1]
                                        fsy[1] = self.sitegeoconl[i][0]
                                        if fsy[0]==fsy[1]:
                                            break;
                                    #find the change in energy due to the reaction
                                    #1.) Atomization - same as before
                                    for il in range(len(ist)):
                                        sind = 0;
                                        if ist[il][0]=='g':
                                           sind = 1;
                                        cmpp = ist[il][sind:]
                                        fndo = False
                                        for fl in range(len(fst)):
                                            sind = 0;
                                            if fst[fl][0]=='g':
                                               sind = 1;
                                            cmpf = fst[fl][sind:]
                                            if cmpf==cmpp:
                                                fndo = True
                                        if fndo==False:
                                            if not(np.isnan(self.prt_info[self.all_prt_dict[cmpp],0])):
                                                delE += self.prt_info[self.all_prt_dict[cmpp],0]
                                    for fl in range(len(fst)):
                                        sind = 0;
                                        if fst[fl][0]=='g':
                                            sind = 1;
                                        cmpp = fst[fl][sind:]
                                        fndo = False
                                        for il in range(len(ist)):
                                            sind = 0;
                                            if ist[il][0]=='g':
                                                sind = 1;
                                            cmpi = ist[il][sind:]
                                            if cmpi==cmpp:
                                                fndo = True
                                        if fndo==False:
                                            if not(np.isnan(self.prt_info[self.all_prt_dict[cmpp],0])):
                                                delE -= self.prt_info[self.all_prt_dict[cmpp],0]
                                    #print sy,fsy,delE
                                    #raw_input()
                                    #2.) Site energies
                                    allowed=True
                                    isyst = np.zeros(2); isyst[:] = -1
                                    fsyst = np.zeros(2); fsyst[:] = -1
                                    lput = 0
                                    for il in range(len(ist)):
                                        if ist[il][0]!='g':
                                            isyst[lput] = self.site_prt_dict[ist[il]]
                                            if not(np.isnan(self.prt_info[self.all_prt_dict[ist[il]],1+sy[lput]])):
                                                delE -= self.prt_info[self.all_prt_dict[ist[il]],1+sy[lput]]
                                                lput += 1
                                            else:
                                                allowed=False
                                    lput = 0
                                    for fl in range(len(fst)):
                                        if fst[fl][0]!='g':
                                            fsyst[lput] = self.site_prt_dict[fst[fl]]
                                            if not(np.isnan(self.prt_info[self.all_prt_dict[fst[fl]],1+fsy[lput]])):
                                                delE += self.prt_info[self.all_prt_dict[fst[fl]],1+fsy[lput]]
                                                lput += 1
                                            else:
                                                allowed=False
                                    #print "allowed?", allowed
                                    if allowed:
                                        if min(isyst)<0 or min(fsyst)<0:
                                            print "Error: initial or final state not found"
                                            exit(1)
                                        aeng = float(ae[0])
                                        from sympy.functions.special.delta_functions import Heaviside
                                        delEf= Heaviside(delE)*abs(delE)+aeng
                                        delEb= Heaviside(-delE)*abs(delE)+aeng
                                        rf   = lambda T, dE=delEf, pe=pf: pe(T)*np.exp(-float(dE/(T*R*0.001)))
                                        rb   = lambda T, dE=delEb, pe=pb: pe(T)*np.exp(-float(dE/(T*R*0.001)))
                                        #print sy, fsy, delE, aeng, delEf, delEb, rf, rb, pf, pb
                                        #print isyst, fsyst
                                        #raw_input()
                                        #Order the reaction
                                        import copy
                                        oisy = copy.deepcopy(sy);
                                        ofsy = copy.deepcopy(fsy);
                                        if oisy[0]>oisy[1]:
                                            tmp=oisy[0];  oisy[0]=oisy[1];   oisy[1]=tmp;
                                            tmp=isyst[0]; isyst[0]=isyst[1]; isyst[1]=tmp;
                                        if ofsy[0]!=oisy[0]:
                                            tmp=ofsy[0];   ofsy[0]=ofsy[1];   ofsy[1]=tmp;
                                            tmp=fsyst[0];  fsyst[0]=fsyst[1]; fsyst[1]=tmp;
                                        if oisy[0]==oisy[1]:
                                            if isyst[0]>isyst[1]:
                                                tmp=isyst[0]; isyst[0]=isyst[1]; isyst[1]=tmp;
                                            if fsyst[0]>fsyst[1]:
                                                tmp=fsyst[0];  fsyst[0]=fsyst[1]; fsyst[1]=tmp;
                                        #print sy, fsy, delE, aeng, delEf, delEb, rf, rb, pf, pb
                                        #print isyst, fsyst
                                        pairind = self.sitegeocon[(oisy[0],oisy[1])]
                                        if(not(oisy[0]!=oisy[1] and isyst[0]==isyst[1] and fsyst[0]==fsyst[1] and (j!=0 or k!=0)) and
                                           not(oisy[0]!=oisy[1] and isyst[0]==isyst[1] and j!=0) and
                                           not(oisy[0]!=oisy[1] and fsyst[0]==fsyst[1] and k!=0)):
                                            if(oisy[0]!=oisy[1]):
                                                scl=0.5
                                                if(fsyst[0]==fsyst[1]):
                                                    if(rf(self.Temp)>cu):
                                                        pairrx.append([pairind,isyst[0],isyst[1],fsyst[0],fsyst[1],rf])
                                                else:
                                                    rfs = lambda T, rr=rf, sc=scl: rr(T)*sc
                                                    if(rf(self.Temp)>cu):
                                                        pairrx.append([pairind,isyst[0],isyst[1],fsyst[0],fsyst[1],rfs])
                                                if(isyst[0]==isyst[1]):
                                                    if(rb(self.Temp)>cu):
                                                        pairrx.append([pairind,fsyst[0],fsyst[1],isyst[0],isyst[1],rb])
                                                else:
                                                    rbs = lambda T, rr=rb, sc=scl: rr(T)*sc
                                                    if(rb(self.Temp)>cu):
                                                        pairrx.append([pairind,fsyst[0],fsyst[1],isyst[0],isyst[1],rbs])
                                            else: #oisy[0]==oisy[1]
                                                scl=0.5
                                                if(isyst[0]==isyst[1] and fsyst[0]==fsyst[1]):
                                                    if(rf(self.Temp)>cu):
                                                        pairrx.append([pairind,isyst[0],isyst[1],fsyst[0],fsyst[1],rf])
                                                    if(rb(self.Temp)>cu):
                                                        pairrx.append([pairind,fsyst[0],fsyst[1],isyst[0],isyst[1],rb])
                                                elif(isyst[0]!=isyst[1] and fsyst[0]==fsyst[1]):
                                                    if(rf(self.Temp)>cu):
                                                        pairrx.append([pairind,isyst[0],isyst[1],fsyst[0],fsyst[1],rf])
                                                        pairrx.append([pairind,isyst[1],isyst[0],fsyst[0],fsyst[1],rf])
                                                    if(rb(self.Temp)>cu):
                                                        rbs = lambda T, rr=rb, sc=scl: rr(T)*sc
                                                        pairrx.append([pairind,fsyst[0],fsyst[1],isyst[0],isyst[1],rbs])
                                                        pairrx.append([pairind,fsyst[0],fsyst[1],isyst[1],isyst[0],rbs])
                                                elif(isyst[0]==isyst[1] and fsyst[0]!=fsyst[1]):
                                                    if(rf(self.Temp)>cu):
                                                        rfs = lambda T, rr=rf, sc=scl: rr(T)*sc
                                                        pairrx.append([pairind,isyst[0],isyst[1],fsyst[0],fsyst[1],rfs])
                                                        pairrx.append([pairind,isyst[0],isyst[1],fsyst[1],fsyst[0],rfs])
                                                    if(rb(self.Temp)>cu):
                                                        pairrx.append([pairind,fsyst[0],fsyst[1],isyst[0],isyst[1],rb])
                                                        pairrx.append([pairind,fsyst[1],fsyst[0],isyst[0],isyst[1],rb])
                                                else:
                                                    if(rf(self.Temp)>cu):
                                                        rfs = lambda T, rr=rf, sc=scl: rr(T)*sc
                                                        pairrx.append([pairind,isyst[0],isyst[1],fsyst[0],fsyst[1],rfs])
                                                        pairrx.append([pairind,isyst[0],isyst[1],fsyst[1],fsyst[0],rfs])
                                                        pairrx.append([pairind,isyst[1],isyst[0],fsyst[0],fsyst[1],rfs])
                                                        pairrx.append([pairind,isyst[1],isyst[0],fsyst[1],fsyst[0],rfs])
                                                    if(rb(self.Temp)>cu):
                                                        rbs = lambda T, rr=rb, sc=scl: rr(T)*sc
                                                        pairrx.append([pairind,fsyst[0],fsyst[1],isyst[0],isyst[1],rbs])
                                                        pairrx.append([pairind,fsyst[0],fsyst[1],isyst[1],isyst[0],rbs])
                                                        pairrx.append([pairind,fsyst[1],fsyst[0],isyst[0],isyst[1],rbs])
                                                        pairrx.append([pairind,fsyst[1],fsyst[0],isyst[1],isyst[0],rbs])
        siterx.sort()
        pairrx.sort()
        self.sitert  = np.zeros((len(siterx),len(siterx[0])-1), dtype=np.int32)
        self.siterr  = np.zeros( len(siterx),    dtype=float)
        self.siterl  = []
        self.pairrt  = np.zeros((len(pairrx),len(pairrx[0])-1), dtype=np.int32)
        self.pairrr  = np.zeros( len(pairrx),    dtype=float)
        self.pairrl  = []
        for i in range(len(siterx)):
            for j in range(len(siterx[i])-1):
                self.sitert[i,j] = siterx[i][j]
            self.siterr[i] = siterx[i][len(siterx[i])-1](self.Temp)
            self.siterl.append(siterx[i][len(siterx[i])-1])
        for i in range(len(pairrx)):
            for j in range(len(pairrx[i])-1):
                self.pairrt[i,j] = pairrx[i][j]
            self.pairrr[i] = pairrx[i][len(pairrx[i])-1](self.Temp)
            self.pairrl.append(pairrx[i][len(pairrx[i])-1])
        #print self.siterx
        #print self.siterxr
        #print len(siterx), len(pairrx)
        #raw_input()
    def updatesystemTemp(self,Temp):
        for i in range(len(self.siterl)):
            self.siterr[i] = self.siterl[i](Temp)
        for i in range(len(self.pairrl)):
            self.pairrr[i] = self.pairrl[i](Temp)

