#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printmeshcell.c                                */
/*                                                                           */
/* Created:       2015/08/31 (VVa)                                           */
/* Last modified: 2015/08/31 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Prints a mesh cell as shown below, based on the point list   */
/*                                                                           */
/*                                                                           */
/* Comments:    -For debugging of fixhexmesh.c                               */
/*                                                                           */
/*                            3---5                                          */
/*             7-----6        |\ /|                                          */
/*            /|    /|        | 4 |           4            3                 */
/*           4-----5 |        | | |          /|\          /|\                */
/*           | 3---|-2        0-|-2        3-----2       0-|-2               */
/*           |/    |/          \|/        /     /         \|/                */
/*           0-----1            1        0-----1           1                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MeshCellFromCGNS:"

/*****************************************************************************/

void  PrintMeshCell(long hex[8], long celltype) 
{
  
  switch (celltype)
    {
    case MESH_CELL_TYPE_TET:
      fprintf(out,"::\n");
      fprintf(out,"     %2ld \n",hex[3]);
      fprintf(out,"     /|\\     \n");
      fprintf(out,"    / | \\   \n");
      fprintf(out,"   /  |  \\   \n");
      fprintf(out,"  /   |   \\   \n");
      fprintf(out,"%2ld ---|-- %2ld\n",hex[0],hex[2]);
      fprintf(out,"  \\   |   /\n");
      fprintf(out,"   \\  |  /\n");
      fprintf(out,"    \\ | /\n");
      fprintf(out,"     \\|/\n");
      fprintf(out,"     %2ld \n",hex[1]);

      break;

    case MESH_CELL_TYPE_PYRAMID:
      fprintf(out,"::\n");
      fprintf(out,"         %2ld \n",hex[4]);
      fprintf(out,"         /|\\     \n\n");      
      fprintf(out,"    %2ld ------ %2ld\n",hex[3],hex[2]);
      fprintf(out,"    /         /\n");
      fprintf(out,"   /         / \n");
      fprintf(out,"  /         /  \n");
      fprintf(out,"%2ld ------ %2ld\n",hex[0],hex[1]);

      break;

    case MESH_CELL_TYPE_PRISM:
      fprintf(out,"    %2ld ------ %2ld\n",hex[3],hex[5]);
      fprintf(out,"     |\\       /|\n");
      fprintf(out,"     | \\     / |\n");
      fprintf(out,"     |  \\   /  |\n");
      fprintf(out,"     |   %2ld    |\n",hex[4]);
      fprintf(out,"     |    |    |\n");      
      fprintf(out,"     |    |    |\n");
      fprintf(out,"    %2ld----|---%2ld\n",hex[0],hex[2]);
      fprintf(out,"      \\   |   /\n");
      fprintf(out,"       \\  |  / \n");
      fprintf(out,"        \\ | /  \n");
      fprintf(out,"         \\|/   \n");
      fprintf(out,"          %2ld  \n",hex[1]);

      break;

    case MESH_CELL_TYPE_HEX:
      fprintf(out,"    %2ld ------ %2ld\n",hex[7],hex[6]);
      fprintf(out,"    /|        /|\n");
      fprintf(out,"   / |       / |\n");
      fprintf(out,"  /  |      /  |\n");
      fprintf(out,"%2ld ------ %2ld   |\n",hex[4],hex[5]);
      fprintf(out," |   |     |   |\n");      
      fprintf(out," |  %2ld ----|- %2ld\n",hex[3],hex[2]);
      fprintf(out," |  /      |  /\n");
      fprintf(out," | /       | / \n");
      fprintf(out," |/        |/  \n");
      fprintf(out,"%2ld ------ %2ld\n",hex[0],hex[1]);

      break;
    default:

      Die(FUNCTION_NAME, "Invalid cell type %ld", celltype);

    }

}
#ifdef __cplusplus 
} 
#endif 
