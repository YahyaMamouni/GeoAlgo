#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <stack>


#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

typedef Polyhedron::Facet_const_iterator Facet_iterator;
typedef Polyhedron::Vertex_const_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_const_iterator Halfedge_iterator;
typedef Polyhedron::Halfedge_around_facet_const_circulator Halfedge_facet_circulator;

typedef std::map<Polyhedron::Facet_const_handle, double> Facet_double_map;
typedef std::map<Polyhedron::Facet_const_handle, int> Facet_int_map;

double table[10][3] = {
    {1.0, 0.0, 0.0},    // red
    {0.0, 1.0, 0.0},    // green
    {0.0, 0.0, 1.0},    // blue
    {1.0, 1.0, 0.0},    // yellow
    {1.0, 0.0, 1.0},    // magenta
    {0.0, 1.0, 1.0},    // cyan
    {0.5, 0.5, 0.5},    // gray
    {1.0, 0.5, 0.0},    // orange
    {0.5, 0.0, 1.0},    // purple
    {0.0, 0.5, 1.0}     // sky blue
};

/// @brief map all the values from [min, max] to [0, 1]
/// @param facetMap non-const reference to the map (it is an in/out parameter)
void normalizeMap(Facet_double_map &facetMap)
{
	double maxValue = facetMap.begin()->second;
	double minValue = facetMap.begin()->second;

	// look for min and max value in the map
	for (const auto &elem : facetMap)
	{
		if (elem.second > maxValue)
		{
			maxValue = elem.second;
		}
		if (elem.second < minValue)
		{
			minValue = elem.second;
		}
	}

	for (auto &elem : facetMap)
	{
		elem.second -= minValue;
		elem.second /= (maxValue-minValue);
	}
}

/// @brief Generate in a .off file a colored mesh according to a value map (green to red shades)
/// @param mesh the input mesh
/// @param facetMap map of values between 0 and 1 (see "normalize()") for each facet of mesh
/// @param filePath path to the colored .off file to be generated
void writeOFFfromValueMap(const Polyhedron& mesh, const Facet_double_map& facetMap, std::string filePath)
{
	std::ofstream in_myfile;
	in_myfile.open(filePath);

	CGAL::set_ascii_mode(in_myfile);

	in_myfile << "COFF" << std::endl // "COFF" makes the file support color informations
			  << mesh.size_of_vertices() << ' ' 
			  << mesh.size_of_facets() << " 0" << std::endl; 
			  // nb of vertices, faces and edges (the latter is optional, thus 0)

	std::copy(mesh.points_begin(), mesh.points_end(),
			  std::ostream_iterator<Kernel::Point_3>(in_myfile, "\n"));

	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{
		Halfedge_facet_circulator j = i->facet_begin();

		CGAL_assertion(CGAL::circulator_size(j) >= 3);

		in_myfile << CGAL::circulator_size(j) << ' ';
		do
		{
			in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());

		} while (++j != i->facet_begin());

		in_myfile << std::setprecision(5) << std::fixed; //set the format of floats to X.XXXXX

		auto redValue = 1-facetMap.at(i); // low values will be closer to red
		auto greenValue = facetMap.at(i); // high values will be closer to green
		auto blueValue = 0.0;

		in_myfile << " " << redValue << " " << greenValue << " " << blueValue;

		in_myfile << std::endl;
	}

	in_myfile.close();

	std::cout << "Le résultat a été exporté dans " << filePath << " !" << std::endl;
}

/**
 @brief Computes the perimeter of each face in the input mesh and returns a map associating each face with its corresponding perimeter.

 @param mesh The input polyhedral mesh

 @return A map associating each face of the mesh with its corresponding perimeter
*/
Facet_double_map computePerimMap(const Polyhedron &mesh)
{
	Facet_double_map out; // Create an empty map to store the perimeters of the faces
	
	// Loop over each face in the mesh and compute its perimeter
	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{
		double current_perimeter = 0.; // Initialize the current perimeter to zero
		Halfedge_facet_circulator j = i->facet_begin(); // Create a circulator to loop over the halfedges of the current face
		do
		{	// Compute the length of the current halfedge and add it to the current perimeter
			current_perimeter += std::sqrt(CGAL::squared_distance(j->vertex()->point(), j->opposite()->vertex()->point()));
		} while (++j != i->facet_begin());

		std::cout << "perim(" << std::distance(mesh.facets_begin(), i) << ")=" << current_perimeter << std::endl;

		out[i] = current_perimeter;
	}

	return out;
}

Facet_double_map computeAreaMap(const Polyhedron &mesh)
{
	Facet_double_map out;

	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{
	
		Halfedge_facet_circulator j = i->facet_begin();
		
		Polyhedron::Vertex_const_handle firstVertex = j->vertex();

		double current_area = 0;
		// a facet is not necessarily a triangle, so we decompose one facet into multiple triangles,
		// and sum up all their areas. Only works for convex faces.
		// (illustration: http://mathbitsnotebook.com/JuniorMath/Polygons/polygons3g.jpg)
		do
		{
			current_area += CGAL::squared_area(firstVertex->point(), j->vertex()->point(), j->opposite()->vertex()->point()); 
		} while (++j != i->facet_begin());

		std::cout << "area(" << std::distance(mesh.facets_begin(), i) << ")=" << current_area << std::endl;

		out[i] = current_area;
	}

	return out;
}

/**
	@brief Computes the threshold value for a given mesh and values.
	@param mesh The input mesh.
	@param values The input values.
	@return The threshold value.
*/
double seuil(const Polyhedron &mesh,Facet_double_map & values){
	double out;
	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{
		out +=values[i];
	}
	return out/ values.size();
	
}

/**

	Assign each facet of a given mesh to a class based on its value compared to a threshold value.

	@param mesh The mesh to analyze.

	@param values A map of double values associated with each facet of the mesh.

	@return A map of integer values indicating the class of each facet.

	This function computes a threshold value by calling the seuil() function, and then assigns each facet of

	the given mesh to a class based on its value compared to this threshold. Facets with values less than 1/4 of

	the threshold are assigned to class 3, facets with values less than 1/2 of the threshold are assigned to

	class 2, facets with values less than 3/4 of the threshold are assigned to class 1, and facets with values

	greater than or equal to 3/4 of the threshold are assigned to class 0.
*/
Facet_int_map simpleThreshold(Polyhedron& mesh, Facet_double_map& values)
{
    Facet_int_map out;
    double s = seuil(mesh, values); // Compute the threshold value

    std::cout << s << std::endl; // Print the threshold value

    // Assign each facet to a class based on its value compared to the threshold
    for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
    {
        if (values[i] < s / 4) { // If the value is less than 1/4 of the threshold
            out[i] = 3; // Assign it to class 3
        }
        else if (values[i] < s / 2) { // If the value is less than 1/2 of the threshold
            out[i] = 2; // Assign it to class 2
        }
        else if (values[i] < s * 3 / 4) { // If the value is less than 3/4 of the threshold
            out[i] = 1; // Assign it to class 1
        }
        else { // Otherwise (if the value is greater than or equal to 3/4 of the threshold)
            out[i] = 0; // Assign it to class 0
        }
    }

    return out;
}

/**

	@brief Writes an OFF file from the given mesh, assigning each face to a color class based on its integer value in the classes map.

	@param mesh The input mesh.

	@param classes A map that associates each face in the mesh with an integer value representing its color class.

	@param filePath The path to the output file.
*/
void writeOFFfromClasse(const Polyhedron& mesh, const Facet_int_map& classes, std::string filePath)
{
	// Open the output file and set the output mode to ASCII
	std::ofstream in_myfile;
	in_myfile.open(filePath);
	CGAL::set_ascii_mode(in_myfile);

	// Write the header of the OFF file with the number of vertices, faces and edges
	in_myfile << "COFF" << std::endl // "COFF" makes the file support color informations
	<< mesh.size_of_vertices() << ' '
	<< mesh.size_of_facets() << " 0" << std::endl;
	// nb of vertices, faces and edges (the latter is optional, thus 0)

    // Write each facet and its color
	std::copy(mesh.points_begin(), mesh.points_end(),
			  std::ostream_iterator<Kernel::Point_3>(in_myfile, "\n"));

	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{
		Halfedge_facet_circulator j = i->facet_begin();

		CGAL_assertion(CGAL::circulator_size(j) >= 3);

		in_myfile << CGAL::circulator_size(j) << ' ';
		do
		{
			in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());

		} while (++j != i->facet_begin());

		in_myfile << std::setprecision(5) << std::fixed; //set the format of floats to X.XXXXX

		auto redValue = table[classes.at(i)][0]; // low values will be closer to red
		auto greenValue = table[classes.at(i)][1]; // high values will be closer to green
		auto blueValue = table[classes.at(i)][2];

		in_myfile << " " << redValue << " " << greenValue << " " << blueValue;

		in_myfile << std::endl;
	}

	in_myfile.close();

	std::cout << "Le résultat a été exporté dans " << filePath << " !" << std::endl;
}

/*
void marquer_composante_connexe(const Polyhedron& mesh, Polyhedron::Facet_const_handle facet, Facet_int_map& segmentation, int label) 
{   
	segmentation[facet] = label%9; 
	for (Polyhedron::Halfedge_around_facet_const_circulator haf = facet->facet_begin(); haf != haf->facet_begin(); ++haf) 
	{
		Polyhedron::Facet_const_handle neighbor = haf->opposite()->facet(); 
		if (segmentation.find(neighbor) == segmentation.end() && segmentation[facet] == segmentation[neighbor]) 
		{ 
			marquer_composante_connexe(mesh, neighbor, segmentation, label); 
		} 
	} 
} 
// Fonction principale pour produire une nouvelle segmentation avec des identifiants différents pour chaque composante connexe 
Facet_int_map segmentationParCC(Polyhedron & mesh, const Facet_int_map & segmentation) 
{ 	
	Facet_int_map nouvelle_segmentation; int label = 0; 
    for (Facet_iterator f = mesh.facets_begin(); f != mesh.facets_end(); ++f) 
    { 
		if (nouvelle_segmentation.find(f) == nouvelle_segmentation.end()) 
		{ 
			marquer_composante_connexe(mesh, f, nouvelle_segmentation, label); ++label; } 
	} 
		
	return nouvelle_segmentation; 
}*/


/**
	@brief Iterates over all faces in the mesh and updates the segmentation of the faces based on their connectivity.
	@param mesh The Polyhedron mesh to iterate over.
	@param outSegmentation A reference to the Facet_int_map containing the segmentation information for each face in the mesh.
	@param currentFace A Facet_iterator representing the current face being processed.
	@param visitedFaces A boolean array indicating whether each face has been visited during the iteration process.
	@param newClass An integer representing the new class for faces in the mesh that have not been visited before.
	@param previousClass A double representing the class of the current face in the previous segmentation.
*/
void iterateOverFaces(Polyhedron& mesh, Facet_int_map& outSegmentation, Facet_iterator currentFace, bool* visitedFaces, int& newClass, double previousClass) {
    Facet_iterator initialFace = mesh.facets_begin();				// create a reference to the initial face for simplification
    if(!visitedFaces[std::distance(initialFace, currentFace)]) {	// if the current face has not been visited yet
        if(outSegmentation[currentFace] == previousClass) {			// if its class in the new segmentation is the same as in the previous one
            visitedFaces[std::distance(initialFace, currentFace)] = true;	// indicate that the face has been visited
            outSegmentation[currentFace] = newClass;						// update the class in the segmentation table
            std::cout << "class(" << std::distance(initialFace, currentFace) << ")=" << outSegmentation[currentFace] << std::endl;

            Halfedge_facet_circulator neighborFace = currentFace->facet_begin();
            do {
                iterateOverFaces(mesh, outSegmentation, neighborFace->opposite()->facet(), visitedFaces, newClass, previousClass); // do the same thing for the neighboring class
            } while(++neighborFace != currentFace->facet_begin());
        }
    }
}


/**
	@brief Segments the mesh into connected components using an iterative algorithm.
	This function takes a mesh and a facet segmentation as input and segments the mesh into connected components
	using an iterative algorithm. It returns a new facet segmentation where all facets belonging to the same connected
	component are assigned the same integer value.
	@param mesh The input mesh to segment.

	@param segmentation The initial facet segmentation.

	@return The new facet segmentation.
*/
Facet_int_map segmentationParCC(Polyhedron & mesh, Facet_int_map segmentation){

    Facet_int_map output = segmentation; 							                // copy segmentation to modify it
    bool visited[mesh.size_of_facets()] = {false}; 					                // initialize a bool array to perform the traversal
    Facet_iterator initial_face = mesh.facets_begin();					            // create a reference to the initial face to simplify
    double current_segment = segmentation[initial_face]; 				            // segment of the initial face according to the previous segmentation
    int new_segment = 0;
    for(Facet_iterator f = mesh.facets_begin(); f != mesh.facets_end(); ++f){ 	    // iterate over each face
        
        if(!visited[std::distance(initial_face,f)]){ 						        // if it has not been visited
            new_segment++;														    // create a new segment
            current_segment = segmentation[f]; 							            // keep the value of its segment in the previous segmentation
            iterateOverFaces(mesh, output, f, visited, new_segment, current_segment);	// launch the iterative traversal
        }
    }
    return output;
}




int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		std::cerr << "Il manque un paramètre au programme. Veuillez lui donner en entrée un nom de fichier au format off." << std::endl;
		return 1;
	}

	Polyhedron mesh;

	std::ifstream input(argv[1]);

	if (!input || !(input >> mesh) || mesh.is_empty())
	{
		std::cerr << "Le fichier donné n'est pas un fichier off valide." << std::endl;
		return 1;
	}

	auto mapPerim = computePerimMap(mesh);

	auto mapArea = computeAreaMap(mesh);

	auto mapClasses = simpleThreshold(mesh, mapArea);

	auto newClasses = segmentationParCC(mesh,mapClasses); 

	writeOFFfromClasse(mesh, newClasses, argc>=3?argv[2]:"result.off");

	//normalizeMap(mapPerim);

	//writeOFFfromValueMap(mesh, mapPerim, argc>=3?argv[2]:"result.off");

	
	return 0;
}
