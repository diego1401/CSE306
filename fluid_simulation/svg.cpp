#include <vector>
#include "../Raytracer_Assignement1/classes/vector.cpp"

// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name
class Polygon {  
public:
	std::vector<Vector> vertices;
	bool is_wet;
	Polygon(){
		this->is_wet = false;
	};
	Polygon(std::vector<Vector> _vertices){
		this->vertices = _vertices;
		this->is_wet = false;
	}

	double Area2D(){
		double res = 0;
		Vector curr,nxt;
		for(int i=0;i<this->vertices.size();i++){
			curr = this->vertices[(i>0)?(i-1):this->vertices.size()-1];
            nxt  = this->vertices[i]; 
			
			res += curr[0]*nxt[1] - nxt[0]*curr[1];
		}
		return 0.5 * std::abs(res);
	}

	Vector Centroid2d(){
		double A = this->Area2D();
		double Cx = 0; double Cy = 0;
		Vector curr,nxt;
		for(int i=0;i<this->vertices.size();i++){
			curr = this->vertices[(i>0)?(i-1):this->vertices.size()-1];
            nxt = this->vertices[i]; 
			Cx += (curr[0] + nxt[0]) * (curr[0]*nxt[1] - nxt[0]*curr[1]);
			Cy += (curr[1] + nxt[1]) * (curr[0]*nxt[1] - nxt[0]*curr[1]);
		}
		return Vector(Cx,Cy,0)* (1/(6*A));
	}
};	

// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &polygons, std::string filename,Polygon polygon, std::string fillcol) {
		FILE* f = fopen(filename.c_str(), "w+"); 
		fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
		for (int j = 0; j < polygon.vertices.size(); j++) {
			if(polygon.vertices[j].ret_is_liquid()){
			fprintf(f, "<circle cx=\"%3.3f\" cy=\"%3.3f\" r=\"1\" stroke=\"black\" stroke-width=\"1\" fill = \"red\" />", (polygon.vertices[j][0]* 1000), (1000 - polygon.vertices[j][1] * 1000));
			}
			else{
            fprintf(f, "<circle cx=\"%3.3f\" cy=\"%3.3f\" r=\"1\" stroke=\"black\" stroke-width=\"1\" fill = \"black\" />", (polygon.vertices[j][0]* 1000), (1000 - polygon.vertices[j][1] * 1000));
			}

        	}
        for (int i=0; i<polygons.size(); i++) {
		    fprintf(f, "<g>\n");
		    fprintf(f, "<polygon points = \""); 
		    for (int j = 0; j < polygons[i].vertices.size(); j++) {
			    fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
		    }
			// std::cout << polygons[i].is_wet << std::endl;
			if(polygons[i].is_wet){
				fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
		    	fprintf(f, "</g>\n");
			}
			else{
		    fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", "none");
		    fprintf(f, "</g>\n");
			}
        }
		fprintf(f, "</svg>\n");
		fclose(f);
	}


// Adds one frame of an animated svg file. frameid is the frame number (between 0 and nbframes-1).
// polygons is a list of polygons, describing the current frame.
// The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg_animated(const std::vector<Polygon> &polygons, std::string filename, int frameid, int nbframes) {
		FILE* f;
		if (frameid == 0) {
			f = fopen(filename.c_str(), "w+");
			fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
			fprintf(f, "<g>\n");
		} else {
			f = fopen(filename.c_str(), "a+");
		}
		fprintf(f, "<g>\n");
		for (int i = 0; i < polygons.size(); i++) {
			fprintf(f, "<polygon points = \""); 
			for (int j = 0; j < polygons[i].vertices.size(); j++) {
				fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000-polygons[i].vertices[j][1] * 1000));
			}
			fprintf(f, "\"\nfill = \"none\" stroke = \"black\"/>\n");
		}
		fprintf(f, "<animate\n");
		fprintf(f, "	id = \"frame%u\"\n", frameid);
		fprintf(f, "	attributeName = \"display\"\n");
		fprintf(f, "	values = \"");
		for (int j = 0; j < nbframes; j++) {
			if (frameid == j) {
				fprintf(f, "inline");
			} else {
				fprintf(f, "none");
			}
			fprintf(f, ";");
		}
		fprintf(f, "none\"\n	keyTimes = \"");
		for (int j = 0; j < nbframes; j++) {
			fprintf(f, "%2.3f", j / (double)(nbframes));
			fprintf(f, ";");
		}
		fprintf(f, "1\"\n	dur = \"5s\"\n");
		fprintf(f, "	begin = \"0s\"\n");
		fprintf(f, "	repeatCount = \"indefinite\"/>\n");
		fprintf(f, "</g>\n");
		if (frameid == nbframes - 1) {
			fprintf(f, "</g>\n");
			fprintf(f, "</svg>\n");
		}
		fclose(f);
	}