#include <stdio.h>
#include <stdlib.h>

int find_ind(int lengmt, int* part_indxs, int locstate, int istate, int nst, int st_glind)
{
  int ind = 0;
  for(int l=0; l<lengmt; l++)
    if(l!=locstate)
      ind = ind*nst+part_indxs[l];
    else
      ind = ind*nst+istate;
  return ind+st_glind;
}

int find_ind2(int lengmt, int* part_indxs, int locstate1, int istate1, int locstate2, int istate2, int nst, int st_glind){
  int ind = 0;
  for(int l=0; l<lengmt; l++)
    if(l==locstate1)
      ind = ind*nst+istate1;
    else if(l==locstate2)
      ind = ind*nst+istate2;
    else
      ind = ind*nst+part_indxs[l];
      
  return ind+st_glind;
}

int find_ind_shift_ukl(int ip, int lengmt, int* part_indxs, int nst, int st_glind){
  int ind = 0;
  for(int l=0; l<lengmt; l++)
    if(l==0)
      ind = ind*nst+ip;
    else
      ind = ind*nst+part_indxs[l-1];
  return ind+st_glind;
}

int find_ind_shift_ukl_rp(int ip, int lengmt, int* part_indxs, int nst, int st_glind, int pos, int st){
  int ind = 0;
  for(int l=0; l<lengmt; l++)
    if(l==0)
      ind = ind*nst+ip;
    else if(l==pos)
      ind = ind*nst+st;
    else
      ind = ind*nst+part_indxs[l-1];
  return ind+st_glind;
}


int find_ind_shift_ukr(int ip, int lengmt, int* part_indxs, int nst, int st_glind){
  int ind = 0;
  for(int l=0; l<lengmt; l++)
    if(l==lengmt-1)
      ind = ind*nst+ip;
    else
      ind = ind*nst+part_indxs[l+1];
  return ind+st_glind;
}

int find_ind_shift_ukr_rp(int ip, int lengmt, int* part_indxs, int nst, int st_glind, int pos, int st)
{
  int ind = 0;
  for(int l=0; l<lengmt; l++)
    if(l==lengmt-1)
      ind = ind*nst+ip;
    else if(l==pos)
      ind = ind*nst+st;
    else
      ind = ind*nst+part_indxs[l+1];
  return ind+st_glind;
}

void recurse_dmt(int ind_rec, int lengmt, int* part_indxs, int* sytys, double* mt, double* dmt, int nst, 
                 int* sitert, double* siterr, int srrl, int* pairrt, double* pairrr, int prrl, int* connect,
                 double * pairp)
{
  for(int k = 0; k<nst; k++){
    part_indxs[ind_rec]=k;
    sytys[ind_rec]     =1; // cus - hack for now
    if(ind_rec<lengmt-1){
      recurse_dmt(ind_rec+1,lengmt,part_indxs, sytys, mt,dmt,nst, sitert, siterr, srrl, pairrt, pairrr, prrl, connect, pairp);
    }
    else{
      int gl_ind = 0;
      for(int l=0; l<lengmt; l++)
        gl_ind = gl_ind*3+part_indxs[l];
      gl_ind += 18;
      dmt[gl_ind] = 0;

      int i = 0; //reaction index
      int sy= 1; //ALL THE SITE TYPES ARE CUS TYPE FOR NOW in the mixed tiling- simple 'hacky' code
      int ind;

      //site reactions
      while(sitert[3*i]  <  sy && i<srrl){ i++; }
      while(sitert[3*i]  <sy+1 && i<srrl){
        for(int pi=0; pi<lengmt; pi++){//loop through part_indxs
          if(sitert[3*i+1]==part_indxs[pi])
            dmt[gl_ind]-=mt[gl_ind]*siterr[i];
          else if(sitert[3*i+2]==part_indxs[pi]){
            ind = find_ind(lengmt, part_indxs, pi, sitert[3*i+1], nst, 18);
            dmt[gl_ind]+=mt[ind]*siterr[i];
          }
        }
        i++;
      }
      
      //pair reactions
      i = 0;
      int pairtl;
      int rs1, rs2, ss1, ss2;
      int part1_ind, part2_ind;
      double scl1, scl2, sum;
      while(i<prrl){
        rs1 = connect[2*pairrt[5*i+0]];
        rs2 = connect[2*pairrt[5*i+0]+1];
        //first see if each total pair is effected
        for(int pi=0; pi<lengmt-1; pi++){
          //1 find inner pair type
          ss1 = sytys[pi];
          ss2 = sytys[pi+1];
          pairtl = 2;// in general should be a function: findlocalpairtype(ss1,ss2,connect,npair);
          //
          //2 set local site states
          part1_ind = part_indxs[pi];
          part2_ind = part_indxs[pi+1];
          //
          if(pairrt[5*i+0]==pairtl){
            if(pairrt[5*i+1]==part1_ind && pairrt[5*i+2]==part2_ind)//the reaction is contained within the thread pair (initial cond)
              dmt[gl_ind]-=mt[gl_ind]*pairrr[i];
            else if(pairrt[5*i+3]==part1_ind && pairrt[5*i+4]==part2_ind){//the reaction is contained within the thread pair (final cond)
              ind = find_ind2(lengmt, part_indxs, pi, pairrt[5*i+1], pi+1, pairrt[5*i+2], nst, 18);
              dmt[gl_ind]+=mt[ind]*pairrr[i];
            }
          }
          //next accound for influences of mixed tiles
          //hacky, but will work for now:
          if(pairrt[5*i+0]==2){
            if(pi==0)         { scl1 = 1; }
            else              { scl1 = 0; }
            if(pi+1==lengmt-1){ scl2 = 1; }
            else              { scl2 = 0; }
          }                                               //cc
          else if(pairrt[5*i+0]==1){ scl1 = 2; scl2 = 2; }//bc
          else                     { scl1 = 0; scl2 = 0; }//bb
          // now go through the cus level interactions
          if(pairrt[5*i+0]==2){
            if(pi==0){
              sum=0;
              for(int ip=0; ip<nst; ip++)
                sum+=mt[find_ind_shift_ukl_rp(ip,lengmt,part_indxs,nst, 18, 1, pairrt[5*i+1])];

              if(sum>0.0000000001){
                if(pairrt[5*i+1]==part1_ind)
                  dmt[gl_ind]-=scl1*mt[gl_ind]*mt[find_ind_shift_ukl(pairrt[5*i+2],lengmt,part_indxs,nst, 18)]*pairrr[i]/sum;
                else if(pairrt[5*i+3]==part1_ind){
                  ind = find_ind(lengmt, part_indxs, 0, pairrt[5*i+1], nst, 18);
                  dmt[gl_ind]+=scl1*mt[ind]*mt[find_ind_shift_ukl_rp(pairrt[5*i+2],lengmt,part_indxs,nst, 18, 1, pairrt[5*i+1])]*pairrr[i]/sum;
                }
              }
            }if(pi==lengmt-2){
              sum=0;
              for(int ip=0; ip<nst; ip++)
                sum+=mt[find_ind_shift_ukr_rp(ip,lengmt,part_indxs,nst, 18, lengmt-2, pairrt[5*i+1])];

              if(sum>0.0000000001){
                if(pairrt[5*i+1]==part2_ind)
                  dmt[gl_ind]-=scl2*mt[gl_ind]*mt[find_ind_shift_ukr(pairrt[5*i+2],lengmt,part_indxs,nst, 18)]*pairrr[i]/sum;
                else if(pairrt[5*i+3]==part2_ind){
                  ind = find_ind(lengmt, part_indxs, lengmt-1, pairrt[5*i+1], nst, 18);
                  dmt[gl_ind]+=scl2*mt[ind]*mt[find_ind_shift_ukr_rp(pairrt[5*i+2],lengmt,part_indxs,nst, 18, lengmt-2, pairrt[5*i+1])]*pairrr[i]/sum;
                }
              }
            }
          }  
          //bc side lengths
          else if(pairrt[5*i+0]==1){
            sum=0;
            for(int ip=0; ip<nst; ip++)
              sum+=pairp[pairrt[5*i+2]+nst*(ip+nst*pairrt[5*i+0])];//sum over all bridge sites - hacky, but will be fine for now
            if(sum>0.0000000001){
              if(pairrt[5*i+2]==part1_ind)//if the initial cus site matches, then take away
                dmt[gl_ind]-=scl1*mt[gl_ind]*pairrr[i]*pairp[pairrt[5*i+2]+nst*(pairrt[5*i+1]+nst*pairrt[5*i+0])]/sum;
              else if(pairrt[5*i+4]==part1_ind){
                ind = find_ind(lengmt, part_indxs, pi, pairrt[5*i+2], nst, 18);
                dmt[gl_ind]+=scl1*mt[ ind ] *pairrr[i]*pairp[pairrt[5*i+2]+nst*(pairrt[5*i+1]+nst*pairrt[5*i+0])]/sum;
              }
              //repeat for second site if I am at the last index
              if(pi==lengmt-2){
                if(pairrt[5*i+2]==part2_ind)//if the initial cus site matches, then take away
                  dmt[gl_ind]-=scl1*mt[gl_ind]*pairrr[i]*pairp[pairrt[5*i+2]+nst*(pairrt[5*i+1]+nst*pairrt[5*i+0])]/sum;
                else if(pairrt[5*i+4]==part2_ind){
                  ind = find_ind(lengmt, part_indxs, pi+1, pairrt[5*i+2], nst, 18);
                  dmt[gl_ind]+=scl1*mt[ ind ] *pairrr[i]*pairp[pairrt[5*i+2]+nst*(pairrt[5*i+1]+nst*pairrt[5*i+0])]/sum;
                }
              }
            }
          }
        }
        i++;
      }
    }
  }
}

extern "C" void calc_dmt(double * dmt,   // f(y),where y'=f(y)
                         double * mt,    // y,   where y'=f(y)
                         double * dpairp,   // f(y),where y'=f(y)
                         double * pairp,    // y,   where y'=f(y)
                         int    * sitert,   // site reaction types
                         double * siterr,   // site reaction rates
                         int    * pairrt,   // pair reaction types
                         double * pairrr,   // pair reaction rates
                         int    * nopairs,  // no of sites in lat
                         double * pairmtx,  // no pairs per
                         int    * connect,  // 2d array of pair con
                         int lengmt,        // length of the mixed tile
                         int srrl,int prrl, // no of s/p rx types
                         int npty, int nst)  // no of pair types and site states
  {
    int pairty_i;
    int part1_ind;
    int part2_ind;
    int no_pair   = npty;
    int no_part   = nst;
    int gl_ind;
    int sy1;       
    int sy2;   

    int   i;
    int   ind,rs1,rs2;
    double sum,scl1,scl2;
    double nopair,nopair_o1,nopair_o2;

    //transfer mt to pairp
    for(int j=0;  j<18; j++)//replace 27 with 18
      pairp[j] = mt[j];
    //
    int nst_tothe_lengmt_m2 = 1;
    for(int i = 0; i<lengmt-2; i++)
      nst_tothe_lengmt_m2 *= nst;
    for(int s1=0; s1<3; s1++){
      for(int s2=0; s2<3; s2++){
        pairp[18+3*s2+s1] = 0.0;
        int strind = 18+nst_tothe_lengmt_m2*(s1+s2*nst);
        for(int i = 0; i<nst_tothe_lengmt_m2; i++)
          pairp[18+3*s2+s1] += mt[i+strind];
      }
    }

    //take care of the pairs first
    for(int j=0; j<2; j++){//replace 3 with 2
      pairty_i   = j;
      sy1        = connect[2*pairty_i  ];
      sy2        = connect[2*pairty_i+1];
      for(int k1=0; k1<nst; k1++){
        part1_ind = k1;
        for(int k2=0; k2<nst; k2++){
          i = 0;
          part2_ind = k2;
          gl_ind    = part2_ind+no_part*(part1_ind+no_part*pairty_i);
          dpairp[gl_ind]=0;
          //loop through site reactions - update with match
          while(sitert[3*i]  <sy1   && i<srrl){ i++; }
          while(sitert[3*i]  <sy1+1 && i<srrl){
            if(sitert[3*i+1]==part1_ind)
              dpairp[gl_ind]-=pairp[gl_ind]*siterr[i];
            else if(sitert[3*i+2]==part1_ind){
              ind = part2_ind+no_part*(sitert[3*i+1]+no_part*pairty_i);
              dpairp[gl_ind]+=pairp[ind]*siterr[i];
            }
            i++;
          } i = 0;
          while(sitert[3*i]  <sy2   && i<srrl){ i++; }
          while(sitert[3*i]  <sy2+1 && i<srrl){
            if(sitert[3*i+1]==part2_ind)
              dpairp[gl_ind]-=pairp[gl_ind]*siterr[i];
            else if(sitert[3*i+2]==part2_ind){
              ind = sitert[3*i+1]+no_part*(part1_ind+no_part*pairty_i);
              dpairp[gl_ind]+=pairp[ind]*siterr[i];
            }
            i++;
          }//*/

          //pair reactions
          i = 0;
          while(i<prrl){
            rs1 = connect[2*pairrt[5*i+0]];
            rs2 = connect[2*pairrt[5*i+0]+1];
            //first see if the total pair is effected
            if(pairrt[5*i+0]==pairty_i){
              if(pairrt[5*i+1]==part1_ind && pairrt[5*i+2]==part2_ind)//the reaction is contained within the thread pair (initial cond)
                dpairp[gl_ind]-=pairp[gl_ind]*pairrr[i];
              else if(pairrt[5*i+3]==part1_ind && pairrt[5*i+4]==part2_ind)//the reaction is contained within the thread pair (final cond)
                dpairp[gl_ind]+=pairp[pairrt[5*i+2]+no_part*(pairrt[5*i+1]+no_part*pairrt[5*i+0])]*pairrr[i];
            }//
            //next select (part1,sy1) and find out how the reaction effects this due to uncertain connections
            nopair    = (double) nopairs[pairty_i];
            nopair_o1 = pairmtx[pairrt[5*i+0]        +2*no_pair*pairty_i];
            nopair_o2 = pairmtx[pairrt[5*i+0]+no_pair+2*no_pair*pairty_i];
            scl1      = nopair_o1/nopair;
            scl2      = nopair_o2/nopair;
            if(rs1==sy1){
              sum=0;
              for(int ip=0; ip<no_part; ip++)
                sum+=pairp[ip+no_part*(pairrt[5*i+1]+no_part*pairrt[5*i+0])];
              if(sum>0.0000000001){
                if(pairrt[5*i+1]==part1_ind)
                  dpairp[gl_ind]-=scl1*pairp[gl_ind]*pairp[pairrt[5*i+2]+no_part*(pairrt[5*i+1]+no_part*pairrt[5*i+0])]*pairrr[i]/sum;
                else if(pairrt[5*i+3]==part1_ind){
                  ind = part2_ind+no_part*(pairrt[5*i+1]+no_part*pairty_i);
                  dpairp[gl_ind]+=scl1*pairp[ind   ]*pairp[pairrt[5*i+2]+no_part*(pairrt[5*i+1]+no_part*pairrt[5*i+0])]*pairrr[i]/sum;
                }
              }
            }else if(rs2==sy1){
              sum=0;
              for(int ip=0; ip<no_part; ip++)
                sum+=pairp[pairrt[5*i+2]+no_part*(ip+no_part*pairrt[5*i+0])];
              if(sum>0.0000000001){
                if(pairrt[5*i+2]==part1_ind)
                  dpairp[gl_ind]-=scl1*pairp[gl_ind]*pairp[pairrt[5*i+2]+no_part*(pairrt[5*i+1]+no_part*pairrt[5*i+0])]*pairrr[i]/sum;
                else if(pairrt[5*i+4]==part1_ind){
                  ind = part2_ind+no_part*(pairrt[5*i+2]+no_part*pairty_i);
                  dpairp[gl_ind]+=scl1*pairp[ind   ]*pairp[pairrt[5*i+2]+no_part*(pairrt[5*i+1]+no_part*pairrt[5*i+0])]*pairrr[i]/sum;
                }
              }
            }
            //
            if(rs1==sy2){
              sum=0;
              for(int ip=0; ip<no_part; ip++)
                sum+=pairp[ip+no_part*(pairrt[5*i+1]+no_part*pairrt[5*i+0])];
              if(sum>0.0000000001){
                if(pairrt[5*i+1]==part2_ind)
                  dpairp[gl_ind]-=scl2*pairp[gl_ind]*pairp[pairrt[5*i+2]+no_part*(pairrt[5*i+1]+no_part*pairrt[5*i+0])]*pairrr[i]/sum;
                else 
                if(pairrt[5*i+3]==part2_ind){
                  ind = pairrt[5*i+1]+no_part*(part1_ind+no_part*pairty_i);
                  dpairp[gl_ind]+=scl2*pairp[ind   ]*pairp[pairrt[5*i+2]+no_part*(pairrt[5*i+1]+no_part*pairrt[5*i+0])]*pairrr[i]/sum;
                }
              }
            }else if(rs2==sy2){ 
              sum=0;
              for(int ip=0; ip<no_part; ip++)
                sum+=pairp[pairrt[5*i+2]+no_part*(ip+no_part*pairrt[5*i+0])];
              if(sum>0.0000000001){
                if(pairrt[5*i+2]==part2_ind)
                  dpairp[gl_ind]-=scl2*pairp[gl_ind]*pairp[pairrt[5*i+2]+no_part*(pairrt[5*i+1]+no_part*pairrt[5*i+0])]*pairrr[i]/sum;
                else if(pairrt[5*i+4]==part2_ind){
                  ind = pairrt[5*i+2]+no_part*(part1_ind+no_part*pairty_i);
                  dpairp[gl_ind]+=scl2*pairp[ind   ]*pairp[pairrt[5*i+2]+no_part*(pairrt[5*i+1]+no_part*pairrt[5*i+0])]*pairrr[i]/sum;
                }
              }
            }//
            i++;
          }//*/
        }
      }
    }

    int * part_indxs = (int*)malloc(sizeof(int)*lengmt);
    int * sytys      = (int*)malloc(sizeof(int)*lengmt);
    recurse_dmt(0,lengmt,part_indxs, sytys, mt,dmt,nst, sitert, siterr, srrl, pairrt, pairrr, prrl, connect, pairp);
    free(part_indxs);
    free(sytys);

    //transfer pair data to mt
    for(int j=0;  j<18; j++)//replace 27 with 18
      dmt[j] = dpairp[j];
  }

extern "C" void constrain_mt(double * mt, int bigd){
  for(int i = 0; i<bigd; i++)
    if(mt[i]<0)
      mt[i]=0.0;
    else if(mt[i]>1)
      mt[i]=1.0;
}






