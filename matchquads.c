#include <iostream> 
#include <random> 
#include <memory.h> 
#include <mpi.h> 
#include <math.h> 
#include <assert.h> 
#include <sys/time.h> 
#include <unistd.h> 
#include <fcntl.h>

#define DIFF_THRESH 1e-4

struct quad_set
{
	float *quads;
	size_t sz; //size of array in bytes
	size_t nquads; //number of quads=sz/16
};


int print_quad(float *);
bool is_close_quad(float *, float *);
bool is_close_number(float, float);
bool is_suspect(float *);


/****************************************************
 * load_quads
 * Args: 
 * 	fname: name of the quad file to load
 * 	quad_set: structure to load the quad file into
 *
 * Description:
 * 	The quads created by analysis.py are stored into
 * 	binary files. Each quad is a set of consecutive
 * 	4 byte floating point numbers. Using fread 
 * 	we populate the qs->quads array.
****************************************************/
void load_quads(const char *fname, struct quad_set *qs)
{
	FILE * qfd = fopen(fname, "rb");
	
	fseek(qfd, 0L, SEEK_END);
	qs->sz = ftell(qfd);
	qs->nquads = qs->sz/16;
	rewind(qfd);
	
	qs->quads = (float *) malloc(qs->sz);
	fread( qs->quads, 4, qs->sz, qfd );
	fclose(qfd);
}


/****************************************************
* matchquads
* Args:
* 	qa1: set of quads in test image
* 	qa2: set of quads in reference image. 
* 	rank: the rank of the process
* 	size: the total number of processes
* Description:
*	Element wise comparison of a set of quads as 
*	to determine equivalence within a threshold
* 	(DIFF_TRESH)
 
****************************************************/
void match_quads( struct quad_set *qa1, struct quad_set *qa2, int rank, int size )
{
	int count=0;
	int start = rank*(qa1->nquads/size);
	int stop = start+(qa1->nquads/size);

	for(int idx1=start; idx1<stop; idx1++)
	{
		for(int idx2=0; idx2 < qa2->nquads; idx2++)
		{
			if( is_suspect( &qa1->quads[idx1*4] ) )
				continue;
			if( is_close_quad( &qa1->quads[idx1*4], &qa2->quads[idx2*4] ) )
			{
				printf("%i %i\n", idx1, idx2);
				count++;
			}
		}
	}
	printf("%i\n", count);
}


bool is_suspect( float quad[] )
{
	const float suspect1[] = {0.f, 0.f, 0.f, 0.f};
	const float suspect2[] = {0.f, 0.f, 1.f, 1.f};
	if( is_close_quad(quad, (float *) suspect1) )
		return true;

	else if( is_close_quad(quad, (float *) suspect2) )
		return true;

	return false;
}

bool is_close_quad( float quad1[], float quad2[] )
{
	if ( is_close_number( quad1[0], quad2[0] ) )
		return false;

	if ( is_close_number( quad1[1], quad2[1] ) )
		return false;

	if ( is_close_number( quad1[2], quad2[2] ) )
		return false;

	if ( is_close_number( quad1[3], quad2[3] ) )
		return false;

	return true;
}

bool is_close_number(float num1, float num2)
{
	return (fabs(num1 - num2) > DIFF_THRESH);
}

int print_quad(float *quad)
{
	printf("\n[");
	for(int ii=0; ii<4; ii++)
	{
		printf("%f", quad[ii]);
		if (ii != 3)
		{
			printf(",\t");
		}
	}
	printf("]\n\n");
}

int main(int argc, char ** argv)
{
	struct quad_set qa1, qa2;
	//Read in the quads
	load_quads("pointing0064_merged.bin", &qa1);
	load_quads("skv625064874090.bin", &qa2);

	
	MPI_Init(&argc, &argv);

	int rank, size, ii, jj;
	
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//Do the matching. 
	match_quads(&qa1, &qa2, rank, size);
	
	MPI_Finalize();
	free( qa1.quads );
	free( qa2.quads );
}
