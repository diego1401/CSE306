#include "voronoi_PLE.cpp"

int main(){
    // Vector e1(0.2,0.2,0),e2(0.8,0.2,0),e3(0.1,0.7,0),e4(0.9,0.6,0);
    // Polygon* P = new Polygon;
    // P->vertices = {e1,e2,e4,e3};
    // Vector f1(0.3,0.3,0),f2(0.8,0.3),f3(0.8,0.7),f4(0.3,0.7);
    // Polygon* clip = new Polygon;
    // clip->vertices = {f1,f2,f3,f4};
    // save_svg({polygon_clipping(*P,*clip)},"image_clipped.svg","red");
    // save_svg({*P,*clip},"image.svg");

    Vector e1(0.2,0.2,0),e2(0.8,0.2,0),e3(0.1,0.7,0),e4(0.9,0.6,0),
           e5(0.5,0.5,0),e6(0.4,0.4,0);
    Polygon dataset({e1,e2,e3,e4,e5,e6});
    save_svg({dataset},"polygon.svg");
    Polygon SubjectPolygon = get_bounding_box(dataset);
    save_svg({SubjectPolygon},"image_box.svg");
    std::vector<Polygon> P;
    for(int i=0;i<dataset.vertices.size();i++){
        P.push_back(clip_bissectors(SubjectPolygon,i,dataset));
        break;
    }
    // for(int i=0;i<P.size();i++){
    //     Polygon tmp = P[i];
    //     for(int j=0;j<tmp.vertices.size();j++){
    //         SubjectPolygon.vertices.push_back(tmp.vertices[j]);
    //     }
    // }
    save_svg(P,"image_voronoi.svg");
    return 0;
}