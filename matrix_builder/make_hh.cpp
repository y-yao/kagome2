#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <hps/src/hps.h>
#include <fstream>

unsigned long long site[100],size;
int sites;
std::vector<unsigned long long> bits;

int mz2total(unsigned long long b)
{
int mz2,s;

mz2=-sites;
for(s=0;s<sites;++s)
	if((b&site[s])==0)
		mz2+=2;
	
return mz2;
}


unsigned long long state(unsigned long long b)
{
unsigned long long s1,s2,s;

s1=0;
s2=size-1;
s=(s1+s2)/2;

while(bits[s]!=b)
	{
	if(s==s1)
		return s2;
		
	if(bits[s]>b)
		{
		s2=s;
		s=(s1+s2)/2;
		}
	else
		{
		s1=s;
		s=(s1+s2)/2;
		}
	}
	
return s;
}


int main(int argc,char* argv[])
{
char *edgefile,*hfile;
double mz,zz;
int edges,edge[100][2],e,s,mz2,offdiag,i;
unsigned long long b,s0,s1,offstate[100],row;
FILE *fp;

if(argc==4)
	{
	edgefile=argv[1];
	mz=atof(argv[2]);
	mz2=2*mz;
	hfile=argv[3];
	}
else
	{
	fprintf(stderr,"expected three arguments: edgefile mz hfile\n");
	return 1;
	}
	
fp=fopen(edgefile,"r");

fscanf(fp,"%d%d",&edges,&sites);

for(e=0;e<edges;++e)
	fscanf(fp,"%d%d",&edge[e][0],&edge[e][1]);
	
fclose(fp);

for(s=0;s<sites;++s)
	site[s]=1ULL<<s;

size=0;

for (b = 0; b < 1ULL << sites; ++b) {
  if (mz2total(b)==mz2) {
    bits.push_back(b);
  }
}
size = bits.size();

std::ofstream bits_file("bits.dat", std::ofstream::binary);
hps::to_stream(bits, bits_file);

fp=fopen(hfile,"w");

for(row=0;row<size;++row)
	{
	b=bits[row];
		
	zz=0.25*edges;
	offdiag=0;
	
	for(e=0;e<edges;++e)
		{
		s0=site[edge[e][0]];
		s1=site[edge[e][1]];
		
		if(((b&s0)==0)!=((b&s1)==0))
			{
			zz-=0.5;
			offstate[offdiag++]=state((b^s0)^s1);
			}
		}
        
        int n_nonzero = 1;
        for (i = 0; i < offdiag; i++) {
          if (offstate[i] > row) n_nonzero++;
        }
        	
	fprintf(fp,"%3d   %llu ",n_nonzero,row);
	for(i=0;i<offdiag;++i) {
	  if (offstate[i] > row) fprintf(fp,"%llu ",offstate[i]);
        }
		
	fprintf(fp,"%5.2f\n",zz);
	}
	
fclose(fp);

return 0;
}
			
	
	
	
	
	
