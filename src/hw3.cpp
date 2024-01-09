/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Kyle Gan
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif
#include <vector>
#include <imageIO.h>
#include <glm/glm.hpp>
#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

#define _USE_MATH_DEFINES
#include <math.h>
char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0
unsigned char buffer[HEIGHT][WIDTH][3];

struct Ray
{
    glm::dvec3 o;
    glm::dvec3 d;
    Ray(glm::dvec3& o, glm::dvec3& d) : o(o), d(d) {}
    Ray() { o = glm::dvec3(0, 0, 0); d = glm::dvec3(0, 0, 0); }
};

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];
int areaLightRows = 1;
int areaLightCols = 1;
int areaLights = areaLightCols * areaLightRows;
int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);


std::vector<Ray> primary(double px, double py)
{
    double fovh = std::tan((fov * 0.5) * M_PI / 180);
    int ind = 0;
    std::vector<Ray> ret(4);
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            double x = 2 * (px + 0.25 + 0.5 * i) / (double)WIDTH - 1;
            double y = 2 * (py + 0.25 + 0.5 * j) / (double)HEIGHT - 1;
            double rx = x * fovh * (double)WIDTH / (double)HEIGHT;
            double ry = y * fovh;
            glm::dvec3 d(rx, ry, -1.0);
            ret[ind] = Ray(glm::dvec3(0.0, 0.0, 0.0), glm::normalize(d));
            ind++;
        }
    }
    return ret;
}

std::vector<Light> area(Light L)
{
    std::vector<Light> ret(areaLightCols * areaLightRows);
    double width = 2.0;
    double height = 2.0;
    double deltaX = width / areaLightCols;
    double deltaY = height / areaLightRows;
    int index = 0;
    for (int i = 0; i < areaLightRows; i++)
    {
        for (int j = 0; j < areaLightCols; j++)
        {
            ret[index].position[0] = L.position[0] - (width / 2.0) + (deltaX * j) + (deltaX / 2.0);
            ret[index].position[1] = L.position[1];
            ret[index].position[2] = L.position[2] - (height / 2.0) + (deltaY * i) + (deltaY / 2.0);

            // Set the color and intensity of the point light
            for (int i = 0; i < 3; i++)
            {
                ret[index].color[i] = L.color[i] / areaLights;
            }
            index++;
        }
    }
    return ret;
}


// TAKES IN EMPTY glm::dvec3
bool findIntersection(Ray& ray, Sphere& s, glm::dvec3& p)
{
    glm::dvec3 do2c = (ray.o - glm::dvec3(s.position[0], s.position[1], s.position[2]));
    double b = 2.0 * glm::dot(ray.d, do2c);
    double c = glm::dot(do2c, do2c) - s.radius * s.radius;
    double discriminant = b * b - 4.0 * c;
    double t0 = -1e10;
    double t1 = -1e10;
    double t = -1e10;
    if (discriminant < -(1e-10))  return false;
    double absD = std::abs(discriminant);
    // REFERENCE https://www.realtimerendering.com/raytracing/Ray%20Tracing%20in%20a%20Weekend.pdf
    // BASIC CONCEPT INSPIRATION
    if (absD < 1e-10)
    {
        t0 = (-b + std::sqrt(absD)) * 0.5;
        t1 = (-b - std::sqrt(absD)) * 0.5;
        t = t0;
    }
    else
    {
        // OPTMIZATION USED
        if (b < 0) t0 = (-b + std::sqrt(discriminant)) * (double)0.5;
        else t0 = (-b - std::sqrt(discriminant)) * (double)0.5;
        t1 = c / t0;
        t0 = t0;
        t = t0;
    }
    if (t0 < 0 && t1 < 0) return false;
    if (t0 > t1 && t1 > 0) t = t1;
    else t = t0;
    p = ray.o + ray.d * t;
    return true;
}


bool findIntersection(Ray& ray, Triangle& t, glm::dvec3& p)
{
    glm::dvec3 c0(t.v[0].position[0], t.v[0].position[1], t.v[0].position[2]);
    glm::dvec3 c1(t.v[1].position[0], t.v[1].position[1], t.v[1].position[2]);
    glm::dvec3 c2(t.v[2].position[0], t.v[2].position[1], t.v[2].position[2]);
    glm::dvec3 N = glm::normalize(glm::cross((c1 - c0), (c2 - c0)));
    double prod = glm::dot(N, ray.d);
    if (std::abs(prod) < 1e-10) return false;
    double d = -(glm::dot(N, c0));
    double t0 = -(glm::dot(N, ray.o) + d) / prod;
    if (t0 <= 0) return false;
    p = ray.o + ray.d * t0;
    if (glm::dot(glm::cross((c1 - c0), (p - c0)), N) < 0 || glm::dot(glm::cross((c2 - c1), (p - c1)), N) < 0 || glm::dot(glm::cross((c0 - c2), (p - c2)), N) < 0) return false;
    return true;
}


glm::dvec3 Phong(Sphere& s, Light& l, glm::dvec3& cand)
{
    glm::dvec3 lxy(l.position[0], l.position[1], l.position[2]);
    glm::dvec3 L = glm::normalize(lxy - cand);
    glm::dvec3 N = glm::normalize(cand - glm::dvec3(s.position[0], s.position[1], s.position[2]));
    double prod = glm::dot(L, N);
    if (prod > 1) prod = 1;
    else if (prod < 0) prod = 0;

    glm::dvec3 R = glm::normalize(N * (2.0 * (prod)) - L);
    glm::dvec3 V = glm::normalize(-cand);
    double prod2 = glm::dot(R, V);
    if (prod2 > 1) prod2 = 1;
    else if (prod2 < 0) prod2 = 0;

    glm::dvec3 kd(s.color_diffuse[0], s.color_diffuse[1], s.color_diffuse[2]);
    glm::dvec3 ks(s.color_specular[0], s.color_specular[1], s.color_specular[2]);
    double comp = std::pow(prod2, s.shininess);
    return glm::dvec3(l.color[0] * ((kd.x * prod) + (ks.x * comp)), l.color[1] * ((kd.y * prod) + (ks.y * comp)), l.color[2] * ((kd.z * prod) + (ks.z * comp)));
}

glm::dvec3 Phong(Triangle& t, const Light& l, glm::dvec3& cand)
{
    glm::dvec3 c0(t.v[0].position[0], t.v[0].position[1], t.v[0].position[2]);
    glm::dvec3 c1(t.v[1].position[0], t.v[1].position[1], t.v[1].position[2]);
    glm::dvec3 c2(t.v[2].position[0], t.v[2].position[1], t.v[2].position[2]);
    glm::dvec3 c = cand;

    double c0c1c2 = glm::length(glm::cross((c1 - c0), (c2 - c0)));
    double cc1c2 = glm::length(glm::cross((c1 - c), (c2 - c)));
    double c0cc2 = glm::length(glm::cross((c - c0), (c2 - c0)));

    double a = cc1c2 / c0c1c2;
    double b = c0cc2 / c0c1c2;
    double g = 1.0 - a - b;


    glm::dvec3 N = glm::normalize(glm::dvec3((a * t.v[0].normal[0]) + (b * t.v[1].normal[0]) + (g * t.v[2].normal[0]), (a * t.v[0].normal[1]) + (b * t.v[1].normal[1]) + (g * t.v[2].normal[1]),
        (a * t.v[0].normal[2]) + (b * t.v[1].normal[2]) + (g * t.v[2].normal[2])));


    glm::dvec3 lxy(l.position[0], l.position[1], l.position[2]);
    glm::dvec3 L = glm::normalize(lxy - cand);
    double prod = glm::dot(L, N);
    if (prod > 1) prod = 1;
    else if (prod < 0) prod = 0;

    glm::dvec3 R = glm::normalize(N * (2.0 * (prod)) - L);
    glm::dvec3 V = glm::normalize(-cand);
    double prod2 = glm::dot(R, V);
    if (prod2 > 1) prod2 = 1;
    else if (prod2 < 0) prod2 = 0;

    glm::dvec3 kd((a * t.v[0].color_diffuse[0]) + (b * t.v[1].color_diffuse[0]) + (g * t.v[2].color_diffuse[0]), (a * t.v[0].color_diffuse[1]) + (b * t.v[1].color_diffuse[1]) + (g * t.v[2].color_diffuse[1]),
        (a * t.v[0].color_diffuse[2]) + (b * t.v[1].color_diffuse[2]) + (g * t.v[2].color_diffuse[2]));
    glm::dvec3 ks((a * t.v[0].color_specular[0]) + (b * t.v[1].color_specular[0]) + (g * t.v[2].color_specular[0]), (a * t.v[0].color_specular[1]) + (b * t.v[1].color_specular[1]) + (g * t.v[2].color_specular[1]),
        (a * t.v[0].color_specular[2]) + (b * t.v[1].color_specular[2]) + (g * t.v[2].color_specular[2]));
    double shi = (a * t.v[0].shininess) + (b * t.v[1].shininess) + (g * t.v[2].shininess);
    double comp = std::pow(prod2, shi);
    return glm::dvec3(l.color[0] * ((kd.r * prod) + (ks.r * comp)), l.color[1] * ((kd.g * prod) + (ks.g * comp)), l.color[2] * ((kd.b * prod) + (ks.b * comp)));
}

glm::dvec3 computeSpheres(Ray& r, glm::dvec3& c, glm::dvec3& ix, double& min, int& ind)
{
    glm::dvec3 pixel_c = c;
    ind = -1;
    for (int i = 0; i < num_spheres; i++)
    {
        glm::dvec3 cand = glm::dvec3(0.0, 0.0, (double)-1e10);
        if (findIntersection(r, spheres[i], cand) && (cand.z > min))
        {
            pixel_c = glm::dvec3();
            ind = i;
            for (int j = 0; j < num_lights; j++)
            {
                std::vector<Light> mini = area(lights[j]);
                for (int m = 0; m < areaLights; m++)
                {
                    glm::dvec3 lxy(mini[m].position[0], mini[m].position[1], mini[m].position[2]);
                    Ray sray(cand, glm::normalize(lxy - cand));
                    bool covered = false;
                    for (int s = 0; s < num_spheres; s++)
                    {
                        glm::dvec3 p;
                        if (s != i)
                        {
                            if (findIntersection(sray, spheres[s], p))
                            {
                                if (glm::length(lxy - cand) > glm::length(p - cand))
                                {
                                    covered = true;
                                    break;
                                }
                            }
                        }
                    }
                    for (int t = 0; t < num_triangles; t++)
                    {
                        glm::dvec3 p;
                        if (findIntersection(sray, triangles[t], p))
                        {
                            if (glm::length(lxy - cand) > glm::length(p - cand))
                            {
                                covered = true;
                                break;
                            }
                        }
                    }
                    if (!covered)
                        pixel_c += Phong(spheres[i], mini[m], cand);
                }
            }
            ix = cand;
            min = cand.z;
        }
    }
    return pixel_c;
}

glm::dvec3 computeTriangles(Ray& r, const glm::dvec3& c, glm::dvec3 & ix, double& min, int& ind)
{
    glm::dvec3 col = c;
    ind = -1;
    for (int i = 0; i < num_triangles; i++)
    {
        glm::dvec3 cand = glm::dvec3(0.0, 0.0, -1e10);
        double dist = -1e10;
        if (findIntersection(r, triangles[i], cand) && (cand.z > min))
        {
            col = glm::dvec3();
            ind = i;
            for (int j = 0; j < num_lights; j++)
            {
                std::vector<Light> mini = area(lights[j]);
                for (int m = 0; m < areaLights; m++)
                {
                    glm::dvec3 lxy(mini[m].position[0], mini[m].position[1], mini[m].position[2]);
                    Ray sray(cand, glm::normalize(lxy - cand));

                    bool covered = false;
                    for (int k = 0; k < num_spheres; k++)
                    {
                        glm::dvec3 p;
                        if (findIntersection(sray, spheres[k], p))
                        {
                            if (glm::length(lxy - cand) > glm::length(p - cand))
                            {
                                covered = true;
                                break;
                            }
                        }
                    }
                    for (int t = 0; t < num_triangles; t++)
                    {
                        glm::dvec3 p;
                        if (t != i)
                        {
                            if (findIntersection(sray, triangles[t], p))
                            {
                                if (glm::length(lxy - cand) > glm::length(p - cand))
                                {
                                    covered = true;
                                    break;
                                }
                            }
                        }
                    }
                    if (!covered) col += Phong(triangles[i], mini[m], cand);
                }
            }
            min = cand.z;
            ix = cand;
        }
    }
    return col;
}

glm::dvec3 recursiveComputeColor(Ray& r, int reflects)
{
    if (reflects > 1) return glm::dvec3(0.0, 0.0, 0.0);
    double R = 0.0, G = 0.0, B = 0.0;
    double min = -1e10;
    glm::dvec3 ix;
    glm::dvec3 s = glm::dvec3(1.0, 1.0, 1.0);
    glm::dvec3 ref = glm::dvec3();
    int sphereNo, triangleNo;

    s = computeSpheres(r, s, ix, min, sphereNo);
    s = computeTriangles(r, s, ix, min, triangleNo);

    glm::dvec3 ks;
    Ray reflect;
    if (sphereNo != -1)
    {
        Sphere s = spheres[sphereNo];
        glm::dvec3 N = glm::normalize(ix - glm::dvec3(s.position[0], s.position[1], s.position[2]));
        glm::dvec3 L = -(r.d);
        double prod = glm::dot(L, N);
        if (prod > 1) prod = 1;
        else if (prod < 0) prod = 0;

        glm::dvec3 R = glm::normalize(N * (2.0 * (prod)) - L);
        glm::dvec3 origin = ix + R * (1e-10);
        reflect = Ray(origin, R);
        ks = glm::dvec3(s.color_specular[0], s.color_specular[1], s.color_specular[2]);
    }
    else if (triangleNo != -1)
    {
        Triangle t = triangles[triangleNo];
        glm::dvec3 c0(t.v[0].position[0], t.v[0].position[1], t.v[0].position[2]);
        glm::dvec3 c1(t.v[1].position[0], t.v[1].position[1], t.v[1].position[2]);
        glm::dvec3 c2(t.v[2].position[0], t.v[2].position[1], t.v[2].position[2]);
        glm::dvec3 c = ix;
        glm::dvec3 crossProd;
        double c0c1c2 = glm::length(glm::cross((c1 - c0), (c2 - c0))); 
        double cc1c2 = glm::length(glm::cross((c1 - c), (c2 - c)));
        double c0cc2 = glm::length(glm::cross((c - c0), (c2 - c0)));

        double a = cc1c2 / c0c1c2;
        double b = c0cc2 / c0c1c2;
        double g = 1.0 - a - b;
        glm::dvec3 N = glm::normalize(glm::dvec3((a * t.v[0].normal[0]) + (b * t.v[1].normal[0]) + (g * t.v[2].normal[0]),
            (a * t.v[0].normal[1]) + (b * t.v[1].normal[1]) + (g * t.v[2].normal[1]),
            (a * t.v[0].normal[2]) + (b * t.v[1].normal[2]) + (g * t.v[2].normal[2])));
        glm::dvec3 L = -r.d;
        double prod = glm::dot(L, N);
        if (prod > 1) prod = 1;
        else if (prod < 0) prod = 0;

        glm::dvec3 R = glm::normalize(N * (2.0 * (prod)) - L);
        glm::dvec3 origin = ix + (R * (1e-10));
        reflect = Ray(origin, R);
        ks = glm::dvec3((a * t.v[0].color_specular[0]) + (b * t.v[1].color_specular[0]) + (g * t.v[2].color_specular[0]),
            (a * t.v[0].color_specular[1]) + (b * t.v[1].color_specular[1]) + (g * t.v[2].color_specular[1]),
            (a * t.v[0].color_specular[2]) + (b * t.v[1].color_specular[2]) + (g * t.v[2].color_specular[2]));
    }
    if (reflects == 0)
    {
        if (sphereNo == -1 && triangleNo == -1)
            return glm::dvec3(1.0, 1.0, 1.0);
        else
        {
            ref = recursiveComputeColor(reflect, reflects);
            return glm::dvec3((1.0 - 0.5) * s.x + 0.5 * ref.x, (1.0 - 0.5) * s.y + 0.5 * ref.y, (1.0 - 0.5) * s.z + 0.5 * ref.z) + recursiveComputeColor(reflect, reflects + 1) * (0.5 * 2.0);
        }
    }
    else
    {
        if (sphereNo == -1 && triangleNo == -1)
            return glm::dvec3(0.0, 0.0, 0.0);
        else
        {
            ref = recursiveComputeColor(reflect, reflects + 1);
            return glm::dvec3((1.0 - ks.x) * s.x + ks.x * ref.x, (1.0 - ks.y) * s.y + ks.y * ref.y, (1.0 - ks.z) * s.z + ks.z * ref.z) * (1.0 - 0.5) + recursiveComputeColor(reflect, reflects + 2) * 0.5;
        }
    }
}

glm::dvec3 super(int x, int y)
{
    double x0 = 0.0;
    double y0 = 0.0;
    double z0 = 0.0;
    std::vector<Ray> rays = primary(x, y);
    for (int i = 0; i < 4; i++)
    {
        glm::dvec3 c = recursiveComputeColor(rays[i], 0);
        x0 += c.x;
        y0 += c.y;
        z0 += c.z;
    }
    return glm::dvec3(x0/4.0, y0/4.0, x0/4.0);
}


//MODIFY THIS FUNCTION
void draw_scene()
{
  //a simple test output
    for (unsigned int x = 0; x < WIDTH; x++)
    {
        glPointSize(2.0);
        glBegin(GL_POINTS);
        for (unsigned int y = 0; y < HEIGHT; y++)
        {
            glm::dvec3 col = super(x, y);
            col += glm::dvec3(ambient_light[0], ambient_light[1], ambient_light[2]);
            for (int i = 0; i < 3; i++)
            {
                if (col[i] < 0) col[i] = 0;
                else if (col[i] > 1) col[i] = 1;
            }
            plot_pixel(x, y, col.x * 255, col.y * 255, col.z * 255);
        }
        glEnd();
        glFlush();
    }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

