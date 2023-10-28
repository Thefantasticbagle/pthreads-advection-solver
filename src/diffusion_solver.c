/**
 * A solver for the diffusion Problem using pthreads.
 * https://en.wikipedia.org/wiki/Diffusion_equation
 * 
 * @file diffusion_solver.c
 * @author Lars L Ruud
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>

#define _GNU_SOURCE 
#include <pthread.h>

#include "../inc/utils.h"

// Convert 'struct timeval' into seconds in double prec. floating point
#define WALLTIME(t) ((double)(t).tv_sec + 1e-6 * (double)(t).tv_usec)

typedef int64_t int_t;
typedef double real_t;

int counter = 0;
pthread_barrier_t barrier;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t  condition = PTHREAD_COND_INITIALIZER;

int_t
    thread_count = 8,
    y_size,
    x_size,
    iterations,
    snapshot_frequency;

real_t
    *temp[2] = { NULL, NULL },
    *thermal_diffusivity,
    dt;

#define T(x,y)                      temp[0][(y) * (x_size + 2) + (x)]
#define T_next(x,y)                 temp[1][((y) * (x_size + 2) + (x))]
#define THERMAL_DIFFUSIVITY(x,y)    thermal_diffusivity[(y) * (x_size + 2) + (x)]

void* run_simulation ( void* rank );
void time_step ( int_t subgrid_width, int_t subgrid_x0 );
void boundary_condition ( long rank_long, int_t subgrid_width, int_t subgrid_x0 );
void domain_init_global ( void );
void domain_init ( int_t subgrid_width, int_t subgrid_x0 );
void domain_save ( int_t iteration );
void domain_finalize ( void );


void
swap ( real_t** m1, real_t** m2 )
{
    real_t* tmp;
    tmp = *m1;
    *m1 = *m2;
    *m2 = tmp;
}


/**
 * Swaps the values of two arrays.
 * 
 * @param m1 The first array.
 * @param m2 The second array.
*/
int
main ( int argc, char **argv )
{
    // Parse command line arguments
    ARGS *args = parse_args( argc, argv );
    if ( !args )
    {
        fprintf( stderr, "Argument parsing failed\n" );
        exit(1);
    }

    y_size = args->y_size;
    x_size = args->x_size;
    iterations = args->iterations;
    snapshot_frequency = args->snapshot_frequency;

    // Initialize barrier & global varaibles
    pthread_barrier_init ( &barrier, NULL, thread_count );
    domain_init_global();

    // Begin timer
    struct timeval t_start, t_end;
    gettimeofday ( &t_start, NULL );

    // Create threads and start them
    pthread_t* thread_handles = malloc ( thread_count * sizeof(pthread_t) );
    for ( long i = 0; i < thread_count; i++ )
        pthread_create ( &thread_handles[i], NULL, run_simulation, (void*) i );

    // Wait for threads to finish
    for ( long i = 0; i < thread_count; i++ )
        pthread_join ( thread_handles[i], NULL );

    // End timer
    gettimeofday ( &t_end, NULL );
    printf ( "Total elapsed time: %lf seconds\n",
            WALLTIME(t_end) - WALLTIME(t_start)
            );

    // Cleanup
    domain_finalize();
    free ( thread_handles );
    pthread_barrier_destroy ( &barrier );
    pthread_mutex_destroy ( &mutex );
    pthread_cond_destroy ( &condition );
    exit ( EXIT_SUCCESS );
}


/**
 * Waits for all threads to reach this point before continuing.
 * 
 * Note: Standard functions exist for this, but are not universally available.
 * This function is here to show that manual alternatives exist.
*/
void
barrier_manual ( )
{
    pthread_mutex_lock ( &mutex );
    counter++;
    if ( counter >= thread_count )
    {
        counter = 0;
        pthread_cond_broadcast ( &condition );
    } else
        while ( pthread_cond_wait ( &condition, &mutex ) != 0 );
    pthread_mutex_unlock ( &mutex );
}


/**
 * Runs part the simulation as if it was sequential.
 * Used by threads as an entry point.
*/
void*
run_simulation ( void* rank )
{
    long rank_long = (long) rank;

    // Calcualte the thread's boundaries
    int_t   subgrid_width = x_size / thread_count,
            subgrid_x0 = subgrid_width * rank_long;

    if ( rank_long == thread_count - 1 ) 
        subgrid_width += x_size % subgrid_width;

    printf( "rank %li: subrid_width = %li, subgrid_x0 = %li\n", rank_long, subgrid_width, subgrid_x0 );
    
    // Init personal domain
    domain_init ( subgrid_width, subgrid_x0 );

    // Iterate
    for ( int_t iteration = 0; iteration <= iterations; iteration++ )
    {
        // Calculate boundary effects
        // (Data dependency - wait for thread 0 to finish saving and swapping) 
        //barrier_manual();
        pthread_barrier_wait ( &barrier );

        // Simulate a border effects & time-step
        // (No race condition - temp and next are only read or written to, respectively)
        boundary_condition( rank_long, subgrid_width, subgrid_x0 );
        time_step( subgrid_width, subgrid_x0 );

        // (Data dependency - wait for each thread to finish their step before saving/swapping buffers)
        //barrier_manual();
        pthread_barrier_wait ( &barrier );

        // Only one rank needs to do the saving and swapping, as temp/next is shared
        if ( rank_long != 0 )
            continue;

        // Snapshot
        if ( iteration % snapshot_frequency == 0 )
        {
            printf (
                "Iteration %ld of %ld (%.2lf%% complete)\n",
                iteration,
                iterations,
                100.0 * (real_t) iteration / (real_t) iterations
            );

            // Save for next iter
            domain_save ( iteration );
        }

        // Swap buffers for next iter
        swap( &temp[0], &temp[1] );
    }

    // Return
    return NULL;
}


/**
 * Calculates one time-step for a given subgrid.
 * 
 * @param subgrid_width The width of the subgrid.
 * @param subgrid_x0 The leftmost x-coordinate of the grid.
*/
void
time_step ( int_t subgrid_width, int_t subgrid_x0 )
{
    real_t c, t, b, l, r, K, new_value;

    for ( int_t y = 1; y <= y_size; y++ )
    {
        for ( int_t x = subgrid_x0 + 1; x <= subgrid_x0 + subgrid_width; x++ )
        {
            c = T(x, y);

            t = T(x - 1, y);
            b = T(x + 1, y);
            l = T(x, y - 1);
            r = T(x, y + 1);
            K = THERMAL_DIFFUSIVITY(x, y);

            new_value = c + K * dt * ((l - 2 * c + r) + (b - 2 * c + t));

            T_next(x, y) = new_value;
        }
    }
}


/**
 * Handles boundary conditions for a given subgrid.
 * 
 * @param rank_long The rank of the thread responsible for the subgrid.
 * @param subgrid_width The width of the subgrid.
 * @param subgrid_x0 The leftmost x-coordinate of the grid.
*/
void
boundary_condition ( long rank_long, int_t subgrid_width, int_t subgrid_x0 )
{
    for ( int_t x = subgrid_x0 + 1; x <= subgrid_x0 + subgrid_width; x++ )
    {
        T(x, 0) = T(x, 2);
        T(x, y_size+1) = T(x, y_size-1);
    }

    if ( rank_long == 0 )
        for ( int_t y = 1; y <= y_size; y++ )
            T(0, y) = T(2, y);

    else if ( rank_long == thread_count - 1 )
        for ( int_t y = 1; y <= y_size; y++ )
            T(x_size+1, y) = T(x_size-1, y);
}


/**
 * Initializes the global variables of the grid.
 * Remember to finalize the domain when its variables are no longer in use!
 * 
 * @see domain_finalize()
*/
void
domain_init_global ( void )
{
    temp[0] = malloc ( (y_size+2)*(x_size+2) * sizeof(real_t) );
    temp[1] = malloc ( (y_size+2)*(x_size+2) * sizeof(real_t) );
    thermal_diffusivity = malloc ( (y_size+2)*(x_size+2) * sizeof(real_t) );

    dt = 0.1;
}


/**
 * Initializes a thread's part of the domain (a subgrid).
 * Assumes the global grid has been initialized.
 * 
 * @param subgrid_width The width of the subgrid.
 * @param subgrid_x0 The leftmost x-coordinate of the grid.
 * 
 * @see domain_init_global()
*/
void
domain_init ( int_t subgrid_width, int_t subgrid_x0 )
{
    for ( int_t y = 1; y <= y_size; y++ )
    {
        for ( int_t x = subgrid_x0 + 1; x <= subgrid_width + subgrid_x0; x++ )
        {
            real_t temperature = 30 + 30 * sin((x + y) / 20.0);
            real_t diffusivity = 0.05 + (30 + 30 * sin((x_size - x + y) / 20.0)) / 605.0;

            T(x,y) = temperature;
            T_next(x,y) = temperature;
            THERMAL_DIFFUSIVITY(x,y) = diffusivity;
        }
    }
}


/**
 * Saves the grid.
 * Writes to the data/ folder.
 * 
 * @param iteration Which iteration the current state of the grid is for.
*/
void
domain_save ( int_t iteration )
{
    int_t index = iteration / snapshot_frequency;
    char filename[256];
    memset ( filename, 0, 256*sizeof(char) );
    sprintf ( filename, "data/%.5ld.bin", index );

    FILE *out = fopen ( filename, "wb" );
    if ( ! out ) {
        fprintf(stderr, "Failed to open file: %s\n", filename);
        exit(1);
    }

    fwrite( temp[0], sizeof(real_t), (x_size+2)*(y_size+2), out );
    fclose ( out );
}


/**
 * De-allocates the grid from memory.
 * 
 * @see domain_init_global()
*/
void
domain_finalize ( void )
{
    free ( temp[0] );
    free ( temp[1] );
    free ( thermal_diffusivity );
}
