#include <iostream> 
#include <random> 
#include <memory.h> 
#include <mpi.h> 
#include <math.h> 
#include <assert.h> 
#include <sys/time.h> 
#include <unistd.h> 
#include <fcntl.h>

#define DIFF_THRESH 1e-100

struct quad_set
{
	float *quads;
	size_t sz;
};


int print_quad(float *);
bool is_close_quad(float *, float *);
bool is_close_number(float, float);
bool is_suspect(float *);

void load_quads(const char *fname, struct quad_set *qs)
{
	FILE * qfd = fopen(fname, "rb");
	
	fseek(qfd, 0L, SEEK_END);
	qs->sz = ftell(qfd);
	rewind(qfd);
	
	qs->quads = (float *) malloc(qs->sz);
	fread( qs->quads, 4, qs->sz/4, qfd );
	fclose(qfd);
}

void match_quads( struct quad_set *qa1, struct quad_set *qa2 )
{
	int count=0;
	for(int idx1=0; idx1<qa1->sz; idx1=idx1+=4)
	{
		for(int idx2=0; idx2 < qa2->sz; idx2+=4)
		{
			if( is_suspect( &qa1->quads[idx1] ) )
				continue;
			if( is_close_quad( &qa1->quads[idx1], &qa2->quads[idx2] ) )
			{
				printf("%i %i\n", idx1/4, idx2/4);
				count++;
				print_quad( &qa1->quads[idx1] );
				print_quad( &qa2->quads[idx2] );
				return;
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
	return (abs(num1 - num2) > DIFF_THRESH);
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
	load_quads("pointing0064_merged.bin", &qa1);
	load_quads("skv625064874090.bin", &qa2);
	for (int ii=0; ii<24; ii+=4)
		print_quad(&qa1.quads[ii]);
	printf("%li\n", qa1.sz/4);
	//match_quads(&qa1, &qa2);
	free(qa1.quads);
	free(qa2.quads);
}
