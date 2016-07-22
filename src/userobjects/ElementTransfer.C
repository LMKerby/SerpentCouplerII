#include "ElementTransfer.h"
#include "MooseMesh.h"

#include "ElementIntegralVariablePostprocessor.h"

template<>

InputParameters validParams<ElementTransfer>()
{
  InputParameters params = validParams<ElementIntegralVariablePostprocessor>();

  // Define this as an UserObject even though inheriting from Postprocessor
  params.set<std::string>("built_by_action") = "add_user_object";

  return params;
}

// Constructor
ElementTransfer::ElementTransfer(const InputParameters & parameters) :
  ElementIntegralVariablePostprocessor(parameters)
{
}

void ElementTransfer::printNodes() const
{
  int numE = _subproblem.mesh().getMesh().n_active_elem();
  int numN, numS;
  int owner, onEdge;
  std::ifstream in;
  std::ofstream out;
  std::string str;
  std::map <int, int> elemIdMap;
  long totFaces, totOwners, totNeighbors;
  Point nodepoint;
  Node *meshnode;
  const Elem *element;
  const Elem *neighbor;
  std::unique_ptr<Elem> side;
  std::ofstream ifcfile;
  std::ofstream Tfile;
  std::ofstream rhofile;
  std::ofstream mapfile;
  std::ofstream pointsfile;
  std::ofstream facesfile;
  std::ofstream ownersfile;
  std::ofstream neighborsfile;
  std::ofstream facesfile_edge;
  std::ofstream ownersfile_edge;
  std::ofstream *tempfile;

  /* Open the interface file */

  ifcfile.open("mooseifc.in", std::ios::out | std::ios::trunc);

  /* Print the header data */

  ifcfile << "4 fuel1 1\n";
  ifcfile << "mooseifc.out\n";

  numN = _subproblem.mesh().nNodes();

  /* Write number of points (nodes) to interface file */

  ifcfile << numN << std::endl;

  /* Write point (node) coordinates to interface file */

  for(int i = 0 ; i < numN ; i++)\
    {
      meshnode = _subproblem.mesh().nodePtr(i);
      ifcfile << (*meshnode)(0) << " " << (*meshnode)(1) << " " << (*meshnode)(2) << "\n";
    }

  /* Write number of active elements / mesh cells into interface file */

  ifcfile << numE << std::endl;

  int outIdx = 0;
  int tot_num_e = _subproblem.mesh().nElem();
  for(int i = 0 ; i < tot_num_e ; i++)
    {
      /* Get next element */

      element = _subproblem.mesh().elemPtr(i);

      /* Only print active elements! */

      if(!element->active())
	continue;

      /* Get number of nodes on this element */

      numN = element->n_nodes();

      // Get number of faces of this element

      numS = element->n_faces();

      /* Get this elements id (used to get averaged value from array) */

      int eId = element->id();

      /* Print element density, temperature, number of surfaces, output index and cell index */

      ifcfile << "-5.424 " << _average_value_array[eId] <<" " << numS << " " << outIdx << " " << eId << "\n";

      /* Increment output index */

      outIdx++;

      /* Loop over element sides / cell faces */

      for(int j = 0 ; j < numS ; j++)
	{
	  /* Get next side */

	  side = element->side(j);

	  /* Get neighbor element through this side*/

          neighbor = element->neighbor(j);

	  /* Get number of nodes in this side */

	  numN = side->n_nodes();

	  /* Write number of nodes to interface file */

	  ifcfile << numN;

	  /* Write node indices in this side to file */

	  for(int k = 0 ; k < numN ; k++)
	    {
	      ifcfile << " " << side->get_node(k)->id() + 1;
	    }

	  /* Next side will be on next line */

	  ifcfile << std::endl;

	}

    }

  ifcfile.close();

  /***********************************************************/
  /***********************************************************/
  /***********************************************************/
  /* Below this is the print-out to the new interface format */
  /***********************************************************/
  /***********************************************************/
  /***********************************************************/

  /* Open the main interface file */

  ifcfile.open("mooseifc_new.in", std::ios::out | std::ios::trunc);

  /* Print the header data */

  ifcfile << "7 fuel1 1" << std::endl;
  ifcfile << "mooseifc_new.m" << std::endl;
  ifcfile << "1.0 600" << std::endl;
  ifcfile << "5 3 4 4 4" << std::endl;
  ifcfile << "mooseifc_new.points" << std::endl;
  ifcfile << "mooseifc_new.faces" << std::endl;
  ifcfile << "mooseifc_new.owners" << std::endl;
  ifcfile << "mooseifc_new.neighbors" << std::endl;
  ifcfile << "mooseifc_new.rho 1" << std::endl;
  ifcfile << "mooseifc_new.T 1" << std::endl;
  ifcfile << "mooseifc_new.map" << std::endl;

  /* Close the interface file */

  ifcfile.close();

  /****************************************/
  /* Write pointsfile                     */
  /****************************************/

  /* Open the pointsfile */

  pointsfile.open("mooseifc_new.points", std::ios::out | std::ios::trunc);

  /* Get number of points in the mesh */

  numN = _subproblem.mesh().nNodes();

  /* Write number of points (nodes) to points file */

  pointsfile << numN << std::endl;

  /* Write point (node) coordinates to interface file */
  /* With this interface format, Serpent wants the points in meters (not in cm) */

  for(int i = 0 ; i < numN ; i++)\
    {
      meshnode = _subproblem.mesh().nodePtr(i);
      pointsfile << '(' << (*meshnode)(0)/100.0 << " " << (*meshnode)(1)/100.0 << " " << (*meshnode)(2)/100.0 << ")\n";
    }

  /* Close the pointsfile */

  pointsfile.close();

  /*******************************************/
  /* Create the other files at the same time */
  /* when looping over elements              */
  /*******************************************/

  /* The OpenFOAM format requires that the faces that have cells on both sides */
  /* come before faces that have cells on only one side */

  /* Facesfile will contain faces that have a neighbor */
  /* Facesfile_edge will contain neighborless faces (faces at the edge of mesh) */

  facesfile.open("mooseifc_new.faces_noedge", std::ios::out | std::ios::trunc);
  facesfile_edge.open("mooseifc_new.faces_edge", std::ios::out | std::ios::trunc);

  /* Same thing for the owners */

  ownersfile.open("mooseifc_new.owners_noedge", std::ios::out | std::ios::trunc);
  ownersfile_edge.open("mooseifc_new.owners_edge", std::ios::out | std::ios::trunc);

  /* And the neighbors (although there are no neighbors at the edge of the mesh */

  neighborsfile.open("mooseifc_new.neighbors_noedge", std::ios::out | std::ios::trunc);

  /* Temperature, density and output mapping files */
  /* these are easy to write */

  Tfile.open("mooseifc_new.T", std::ios::out | std::ios::trunc);
  rhofile.open("mooseifc_new.rho", std::ios::out | std::ios::trunc);
  mapfile.open("mooseifc_new.map", std::ios::out | std::ios::trunc);

  /* Write number of active elements / mesh cells into some files */

  Tfile << numE << std::endl;
  rhofile << numE << std::endl;
  mapfile << numE << std::endl;

  /* We should print the number of faces first to the file, but we don't know */
  /* that yet so we'll count them */
  /* The situation is the same for the number of owners and neighbors */

  totFaces = 0;
  totOwners = 0;
  totNeighbors = 0;

  /* Reset output index (this is used to identify the elements/cells on Serpent side) */

  outIdx = 0;

  /* Get total number of elements in the mesh */

  tot_num_e = _subproblem.mesh().nElem();

  /* Loop over all elements to print out data */

  for(int i = 0 ; i < tot_num_e ; i++)
    {
      /* Get next element */

      element = _subproblem.mesh().elemPtr(i);

      /* Only print active elements! */

      if(!element->active())
	continue;

      /* Get number of nodes on this element */

      numN = element->n_nodes();

      /* Get number of faces of this element */

      numS = element->n_faces();

      /* Get this elements id (used to get averaged value from array) */

      int eId = element->id();

      /* Map element id to the outIdx, needed when figuring out neighbors */
      /* and owners */

      elemIdMap[eId] = outIdx;

      /* Write temperature to file */

      Tfile << _average_value_array[eId] << std::endl;

      /* Write density to file (could be read from an array) */

      rhofile << "-5.524\n";

      /* Write output mapping index to file */

      mapfile << outIdx + 1 << std::endl;

      /* Increment output index */

      outIdx++;

      /* Loop over element sides / cell faces */

      for(int j = 0 ; j < numS ; j++)
	{
	  /* Get next side */

	  side = element->side(j);

	  /* Get neighbor element through this side*/

          neighbor = element->neighbor(j);

          if (neighbor == NULL)
            {
	      /* This face is at the edge of the geometry */

	      onEdge = 1;

	      /* This element is the owner of the face */

	      owner = 1;

	      /* Increment number of face definitions */

	      totFaces++;
	      totOwners++;
            }
          else
            {
	      /* This face is not at the edge of the geometry */

	      onEdge = 0;

	      /* We'll decide that the owner of the face is the */
	      /* element with the larger id() number */

              if (neighbor->id() < element->id())
		{
		  /* This element is the owner of the face */

		  owner = 1;

		  /* Increment number of face definitions */

		  totFaces++;
		  totOwners++;
		}
              else
		{
		  /* This element is the neighbor of the face */

		  owner = 0;

		  /* Increment number of neighbors */

		  totNeighbors++;
		}
            }

	  /* Write faces out only when processing the owner cell */
	  /* (each face is written out only once) */

	  if (owner)
	    {

	      /********************/
	      /* Write face files */
	      /********************/

	      /* Get number of nodes in this side */

	      numN = side->n_nodes();

	      /* Get file to write to */

	      if (onEdge)
		{
		  tempfile = &facesfile_edge;
		}
	      else
		{
		  tempfile = &facesfile;
		}

	      /* Write number of nodes to interface file */

	      *tempfile << numN << "(";

	      /* Write node indices in this side to file */

	      for(int k = 0 ; k < numN ; k++)
		{
		  *tempfile << " " << side->get_node(k)->id();
		}

	      /* Next side will be on next line */

	      *tempfile << ")" << std::endl;

	      /*********************/
	      /* Write owner files */
	      /*********************/

	      /* Get file to write to */

	      if (onEdge)
		{
		  tempfile = &ownersfile_edge;
		}
	      else
		{
		  tempfile = &ownersfile;
		}

	      /* Write owner index */

	      *tempfile << elemIdMap[element->id()] << "\n";

	      /************************/
	      /* Write neighbor files */
	      /************************/

	      if (!onEdge)
		{
		  /* Write neigbor index */

		  neighborsfile << elemIdMap[neighbor->id()] << "\n";
		}

	    }
	}

    }

  /* Close all of the output files for now */

  Tfile.close();
  rhofile.close();
  mapfile.close();
  facesfile.close();
  facesfile_edge.close();
  ownersfile.close();
  ownersfile_edge.close();
  neighborsfile.close();

  /* We'll still need to do some work with the face, owner and neighbor files */
  /* The total number of instances must be prepended to the final file and the*/
  /* noedge and edge files must be combined*/

  /* Final facesfile */

  out.open("mooseifc_new.faces", std::ios::out | std::ios::trunc);

  /* Write number of faces to faces file */

  out << totFaces << std::endl;

  /* Open file for reading */

  in.open("mooseifc_new.faces_noedge", std::ios::in);

  /* Write non-edge faces */

  if (in.is_open())
  {

    /* Loop over the interface file and get values */

    while(getline(in,str))
      out << str << std::endl;
  }

  in.close();

  /* Open file for reading */

  in.open("mooseifc_new.faces_edge", std::ios::in);

  /* Write edge faces */

  if (in.is_open())
  {

    /* Loop over the interface file and get values */

    while(getline(in,str))
        out<<str;
  }

  in.close();
  out.close();

  /* Final ownersfile */

  out.open("mooseifc_new.owners", std::ios::out | std::ios::trunc);

  /* Write number of owners to owners file */

  out << totOwners << std::endl;

  /* Open file for reading */

  in.open("mooseifc_new.owners_noedge", std::ios::in);

  /* Write non-edge owners */

  if (in.is_open())
  {

    /* Loop over the interface file and get values */

    while(getline(in,str))
        out << str << std::endl;
  }

  in.close();

  /* Open file for reading */

  in.open("mooseifc_new.owners_edge", std::ios::in);

  /* Write edge owners */

  if (in.is_open())
  {

    /* Loop over the interface file and get values */

    while(getline(in,str))
        out << str << std::endl;
  }

  in.close();
  out.close();

  /* Final neighborsfile */

  out.open("mooseifc_new.neighbors", std::ios::out | std::ios::trunc);

  /* Write number of neighbors to neighbors file */

  out << totNeighbors << std::endl;

  /* Open file for reading */

  in.open("mooseifc_new.neighbors_noedge", std::ios::in);

  /* Write non-edge neighbors */

  if (in.is_open())
  {

    /* Loop over the interface file and get values */

    while(getline(in,str))
        out << str << std::endl;
  }

  in.close();
  out.close();

  /* Remove the temporary files */

  std::remove("mooseifc_new.faces_noedge");
  std::remove("mooseifc_new.faces_edge");
  std::remove("mooseifc_new.owners_noedge");
  std::remove("mooseifc_new.owners_edge");
  std::remove("mooseifc_new.neighbors_noedge");

}

void ElementTransfer::initialize()
{
  ElementIntegralVariablePostprocessor::initialize();

  int numE = _subproblem.mesh().nElem();

  /* Initialize the arrays used for calculating the */
  /* element averages */

  _integral_value_array = new Real[numE];
  _volume_value_array = new Real[numE];
  _average_value_array = new Real[numE];

  for(int i = 0 ; i < numE ; i++)
    {
      _integral_value_array[i] = 0.0;
      _volume_value_array[i] = 0.0;
      _average_value_array[i] = 0.0;
    }

}

void ElementTransfer::execute()
{
  // defined in ElementIntegralPostprocessor

  int id = _current_elem->id();

  _integral_value_array[id] += computeIntegral();

  _volume_value_array[id] += _current_elem_volume;

}

void ElementTransfer::threadJoin(const UserObject & y)
{
  ElementIntegralVariablePostprocessor::threadJoin(y);

  /* This should work nicely in parallel */

  const ElementTransfer & otheret = dynamic_cast<const ElementTransfer &>(y);

  int numE = _subproblem.mesh().nElem();

  for(int i = 0 ; i < numE ; i++)
    {
      _integral_value_array[i] += otheret._integral_value_array[i];
      _volume_value_array[i] += otheret._volume_value_array[i];
    }

}

void ElementTransfer::finalize()
{
  int numE = _subproblem.mesh().nElem();

  for(int i = 0 ; i < numE ; i++)
    {
      /* Get the integral of each element */

      gatherSum(_integral_value_array[i]);

      /* Get the volume of each element */

      gatherSum(_volume_value_array[i]);
    }

  for(int i = 0 ; i < numE ; i++)
    {
      /* Average is integral/volume */

      if(_volume_value_array[i] != 0)
	_average_value_array[i] = _integral_value_array[i]/_volume_value_array[i];
      else
	_average_value_array[i] = 0;
    }

}
