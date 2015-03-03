import numpy as np

nst      = 3;
no_ppres = 20;
maxpp    = 50;
minpp    = 0.5;
dlog     = (np.log(maxpp)-np.log(minpp))/no_ppres

def findind(part_indxs,tilesz,stind,nst):
	ind = 0;
	for i in range(tilesz):
		ind = nst*ind+part_indxs[i];
	return ind+stind

def findp(curind,stind,data,len_data,loc,tilesz,nst,state,part_indxs):
	ret = 0;
	if(curind==tilesz):
		ret = data[findind(part_indxs,tilesz,stind,nst)]
	elif(curind!=loc):
		for k in range(nst):
		    part_indxs[curind]=k;
		    ret += findp(curind+1,stind,data,len_data,loc,tilesz,nst,state,part_indxs)
	else:
		part_indxs[curind]=state;
		ret += findp(curind+1,stind,data,len_data,loc,tilesz,nst,state,part_indxs)
	return ret

    

for i in range(2,6):
	foutss     = open(''.join(['results_mixed',str(i),'/sstate.csv']), 'w')
	ss         = np.empty((no_ppres,2+i))# time, bc, and ccc..cc steady states
	part_indxs = np.empty((i), dtype=int)
	for j in range(no_ppres):
		print "Tile size:", i, "pco_indx:", j
		#partial pressure
		pco     = minpp*np.exp(dlog*j)
		#time steps
		fread   = open(''.join(['results_mixed',str(i),'/pairp_pco_',str(j),'.csv']), 'r')
		noline  = 0
		for line in fread:
			noline += 1
		fread.close()
		ProbCO = np.empty((noline,2+i))# there will be 1 result from the bc stuff and i contractions from the cc...cc (i times)

		nst = 3;
		p   = 0;
		fread   = open(''.join(['results_mixed',str(i),'/pairp_pco_',str(j),'.csv']), 'r')
		line_indx  = 0
		for line in fread:
			data = [float(x) for x in line.split()]
			#first calculate P(CO) on bc pairs
			p    =  0
			for ist in range(nst):
				p += data[9+0+nst*ist]
			ProbCO[line_indx,0] = p
			#next calculate P(CO) on each c...c reduction
			for tind in range(i):
				ProbCO[line_indx,1+tind] = findp(0,18,data,len(data),tind,i,nst,0,part_indxs)
			ProbCO[line_indx,i+1] = data[18+nst**i] #record time
			line_indx += 1
		fread.close()
		
		ss[j,0] = pco
		for outind in range(i+1):
			ss[j,1+outind] = ProbCO[noline-1,1+outind] 

		foutd   = open(''.join(['results_mixed',str(i),'/site_pco_',str(j),'.csv']), 'w')
		for outi in range(noline):
			for outj in range(2+i):
				foutd.write(str(ProbCO[outi,outj]))
				foutd.write(" ")
			foutd.write("\n")
		foutd.close()
	#write to ss here
	foutss  = open(''.join(['results_mixed',str(i),'/sstate.csv']), 'w')
	#ss         = np.empty((no_ppres,2+i))# time, bc, and ccc..cc steady states
	for outi in range(no_ppres):
		for outj in range(2+i):
			foutss.write(str(ss[outi,outj]))
			foutss.write(" ")
		foutss.write("\n")
	foutss.close()
	foutss.close()