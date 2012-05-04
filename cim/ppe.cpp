#include "FronTier.h"

int main(int argc, char **argv)
{
	char in_name[100];
	char out_name[100];
	char run_name[100];
	int i;
	FILE *infile,*outfile;
	double L1[10],L2[10],LI[10],L1_order[10],L2_order[10],LI_order[10];

	sprintf(run_name,"%s",argv[1]);
	for (i = 0; i < 4; ++i)
	{
	    sprintf(in_name,"out-%s-%d/run-output",run_name,i+1);
	    infile = fopen(in_name,"r");
	    CursorAfterString(infile,"L1 = ");
	    fscanf(infile,"%lf",&L1[i]);
	    printf("%20.14f\n",L1[i]);
	    CursorAfterString(infile,"L2 = ");
	    fscanf(infile,"%lf",&L2[i]);
	    printf("%20.14f\n",L2[i]);
	    CursorAfterString(infile,"L_inf = ");
	    fscanf(infile,"%lf",&LI[i]);
	    printf("%20.14f\n",LI[i]);
	    if (i != 0)
	    {
		L1_order[i] = log(L1[i-1]/L1[i])/log(2.0);
		L2_order[i] = log(L2[i-1]/L2[i])/log(2.0);
		LI_order[i] = log(LI[i-1]/LI[i])/log(2.0);
	    }
	}
	sprintf(out_name,"%s-error",run_name);
	outfile = fopen(out_name,"w");
	fprintf(outfile,"      L-1       ");
	fprintf(outfile,"      L-2       ");
	fprintf(outfile,"      L-I       ");
	fprintf(outfile,"L1-order ");
	fprintf(outfile," L2-order ");
	fprintf(outfile," LI-order\n");
	for (i = 0; i < 4; ++i)
	{
	    fprintf(outfile,"%14.10f  %14.10f  %14.10f",
			L1[i],L2[i],LI[i]);
	    if (i != 0)
	    	fprintf(outfile,"  %6.4f    %6.4f    %6.4f\n",L1_order[i],
				L2_order[i],LI_order[i]);
	    else
	    	fprintf(outfile,"\n");
	}
}
