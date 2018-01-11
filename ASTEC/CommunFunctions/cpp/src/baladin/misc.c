
#include <misc.h>



void PrintBlock( FILE *f, BLOC* b )
{
  fprintf( f, "(%d,%d,%d), inc=%d, val=%d moy=%f, var=%f diff=%f\n",
	   b->a, b->b, b->c, b->inclus, b->valid, b->moy, b->var, b->diff );
}



void PrintListOfBlocks( FILE *f, BLOCS* b )
{
  int i;
  fprintf( f, "List of blocks: %d valid among %d blocks\n", 
	   b->n_valid_blocks, b->n_blocks );
  for ( i=0; i<b->n_blocks; i++ ) {
    fprintf( f, "#%6d: ", i );
    PrintBlock( f, &(b->bloc[i]) );
  }
}


void PrintListOfValidBlocks( FILE *f, BLOCS* b )
{
  int i, j;
  fprintf( f, "List of blocks: %d valid among %d blocks\n", 
	   b->n_valid_blocks, b->n_blocks );
  for ( i=0, j=0; i<b->n_blocks; i++ ) {
    if ( b->bloc[i].valid == 1 ) {
      fprintf( f, "#%4d(%6d): ", j, i );
      PrintBlock( f, &(b->bloc[i]) );
      j++;
    }
  }
}


void PrintPointersOfBlocks( FILE *f, BLOCS* b )
{
  int i;
  fprintf( f, "Pointers of blocks: %d valid among %d blocks\n", 
	   b->n_valid_blocks, b->n_blocks );
  for ( i=0; i<b->n_blocks; i++ ) {
    fprintf( f, "#%6d: ", i );
    PrintBlock( f, b->p_bloc[i] );
  }
}


void PrintPointersOfValidBlocks( FILE *f, BLOCS* b )
{
  int i, j;
  fprintf( f, "Pointers of blocks: %d valid among %d blocks\n", 
	   b->n_valid_blocks, b->n_blocks );
  for ( i=0, j=0; i<b->n_blocks; i++ ) {
    if ( b->p_bloc[i]->valid == 1 ) {
      fprintf( f, "#%4d(%6d): ", j, i );
      PrintBlock( f, b->p_bloc[i] );
      j++;
    }
  }
}



void PrintParametersOfBlocks( FILE *f, BLOCS* b )
{

  fprintf( f, " valid blocks / blocks = %d\n", b->n_valid_blocks );
  fprintf( f, " allocated blocks = %d\n", b->n_valid_blocks );
  if ( 0 ) {
    fprintf( f, " number of blocks = %d x %d x %d\n", 
	     b->n_blocks_x, b->n_blocks_y, b->n_blocks_z );
    fprintf (f, " block dimensions = %d x %d x %d\n", 
	     b->block_size_x, b->block_size_y, b->block_size_z );
    fprintf (f, " block step = %d x %d x %d\n", 
	     b->block_step_x, b->block_step_y, b->block_step_z );
    fprintf (f, " block border = %d x %d x %d\n", 
	     b->block_border_x, b->block_border_y, b->block_border_z );
  }
}

void PrintParametersForBlocks( FILE *f, PARAM *p )
{
  fprintf ( f, " block dimensions =  %d x %d x %d\n", 
	    p->bl_dx, p->bl_dy, p->bl_dz ); 
  if ( 0 )
    fprintf ( f, " block borders =  %d x %d x %d\n", 
	      p->bl_border_x, p->bl_border_y, p->bl_border_z ); 
  fprintf ( f, " block step  =  %d x %d x %d\n", 
	    p->bl_next_x, p->bl_next_y, p->bl_next_z );

  fprintf ( f, " neighborhood block step =  %d x %d x %d\n", 
	    p->bl_next_neigh_x, p->bl_next_neigh_y, p->bl_next_neigh_z );
  fprintf ( f, " neighborhood dimension =  %d x %d x %d\n", 
	    p->bl_size_neigh_x, p->bl_size_neigh_y, p->bl_size_neigh_z );
  if ( 0 )
    fprintf ( f, "percentage variance = %f (min = %f , dec = %f )\n",
	      p->bl_pourcent_var, p->bl_pourcent_var_min, p->bl_pourcent_var_dec );
  fprintf ( f, "percentage variance = %f (min = %f)\n",
	    p->bl_pourcent_var, p->bl_pourcent_var_min );
}
