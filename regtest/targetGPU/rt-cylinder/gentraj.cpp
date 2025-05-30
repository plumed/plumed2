#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>

double randomn() { return rand()/static_cast<double>(RAND_MAX); }

int main() {

  std::vector<double> tot(8,0), dtot(8,0); bool dens=true;
  unsigned nmols=200, dir; std::vector<double> com(3), dist(3), a1(3), a2(3), len(nmols);
  for(unsigned i=0;i<20;++i){
      if(dens) std::cout<<1+nmols<<std::endl;
      else std::cout<<1+2*nmols<<std::endl;
      std::cout<<10.0<<" "<<10.0<<" "<<10.0<<std::endl;
      std::cout<<"Ar "<<5.0<<" "<<5.0<<" "<<5.0<<std::endl;
      for(unsigned j=0;j<nmols;++j){
         com[0]=10.0*randomn(); com[1]=10.0*randomn(); com[2]=10.0*randomn();
         for(unsigned k=0;k<3;++k) dist[k]=fabs(com[k]-5.0);
         len[j]=randomn(); dir=floor( 3*randomn() );
         if( dist[0]<1.0 ){ tot[1]+=1.0; dtot[1]+=len[j]; }
         if( dist[1]<1.0 ){ tot[2]+=1.0; dtot[2]+=len[j]; }
         if( dist[2]<1.0 ){ tot[3]+=1.0; dtot[3]+=len[j]; }
         if( dist[0]<1.0 && dist[1]<1.0 ){ tot[4]+=1.0; dtot[4]+=len[j]; }
         if( dist[0]<1.0 && dist[2]<1.0 ){ tot[5]+=1.0; dtot[5]+=len[j]; }
         if( dist[1]<1.0 && dist[2]<1.0 ){ tot[6]+=1.0; dtot[6]+=len[j]; }
         if( dist[0]<1.0 && dist[1]<1.0 && dist[2]<1.0 ){ tot[7]+=1.0; dtot[7]+=len[j]; }
         a1[0]=com[0]; a1[1]=com[1]; a1[2]=com[2];
         a2[0]=com[0]; a2[1]=com[1]; a2[2]=com[2];
         a1[dir]-=0.5*len[j]; a2[dir]+=0.5*len[j];
         if(dens){
            std::cout<<"Ar "<<com[0]<<" "<<com[1]<<" "<<com[2]<<std::endl;
         } else {
            std::cout<<"Ar "<<a1[0]<<" "<<a1[1]<<" "<<a1[2]<<std::endl;
            std::cout<<"Ar "<<a2[0]<<" "<<a2[1]<<" "<<a2[2]<<std::endl;
         }
      }
      for(unsigned i=1;i<8;++i){
         if(dens) std::cerr<<"CV "<<i<<" "<<tot[i]<<std::endl; 
         else std::cerr<<"CV "<<i<<dtot[i] / static_cast<double>( nmols )<<std::endl;
         tot[i]=dtot[i]=0.0;
      }
  }

  return 1;
}
