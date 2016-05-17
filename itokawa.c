// -lGLEW
#include <GL/glew.h>
// -lglut
#include <GL/glut.h>
#include <stdio.h>
//-lm
#include <math.h>
// -lgl


struct Trianglew
{
	float p1[4], p2[4], p3[4];
};

void mat_mul(float[16], float[], int, float[]);
void mat_id(float[16]);
void mat_trans(float[16], float, float, float);
void mat_rotx(float[16], float);
void mat_rotx2(float[16], float, float);
void mat_roty(float[16], float);
void mat_roty2(float[16], float, float);
void mat_rotz(float[16], float);
void mat_rotz2(float[16], float, float);
void mat_scale(float[16], float, float, float);
void proyectar(float[16], float[], struct Trianglew);
void dibujar_linea(float a[], int _x1, int _y1, int _x2, int _y2);
void pendiente_menor_menos1(float*, int, int, int, int);
void pendiente_entre_menos1y0(float*, int, int, int, int);
void pendiente_entre_1y0(float*, int, int, int, int);
void pendiente_mayor_1(float*, int, int, int, int);
void calc_mat(float [16]);
void abrir_archivo_y_procesa(struct Trianglew[]);
void setup();
void display();
void draw_pixel(float*, int, int, float, float, float);

int width = 480, height = 480, N=49152;
struct Trianglew triangles[49152];

int main(int argc, char *argv[])
{
	abrir_archivo_y_procesa(triangles);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
	glutInitWindowSize(width,height);
	glutCreateWindow("Itokawa");

	setup();
	glutDisplayFunc(display);
	glutMainLoop();
	return 0;
}

void display(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//glColor3f(0.0f, 0.0f, 0.0f);
	//glRectf(-0.75f,0.75f, 0.75f, -0.75f);
	
	float *a = (float*)malloc(sizeof(float)*height*width*3);

	float m[16];
	calc_mat(m);
	for(int i=0; i<N; ++i){
		proyectar(m, a, triangles[i]);
	}
	glDrawPixels(width, height, GL_RGB, GL_FLOAT, a);
	glutSwapBuffers();
}

void mat_mul(float m1[16], float m2[], int n, float mat[]){
	float column[4];
    for (int j=0; j<n; ++j){
        for (int i=0; i<4; ++i) {
            column[i] = 0;
            for (int l=0; l<4; ++l)
                column[i] += m1[i*4+l] * m2[l*n+j];
        }
        for (int i=0; i<4; ++i)
            mat[i*n+j] = column[i];
    }
    
}

void mat_id(float mat[16]){
    
              mat[1]  = mat[2] = mat[3]  = 
    mat[4]  =           mat[6] = mat[7]  = 
    mat[8]  = mat[9]  =          mat[11] = 
    mat[12] = mat[13] = mat[14]          = 0;
    
    mat[0]  = mat[5]  = mat[10]= mat[15] = 1;

}

void mat_trans(float mat[16], float dx, float dy, float dz){
    float t[16] = {
        1, 0, 0, dx,
        0, 1, 0, dy,
        0, 0, 1, dz,
        0, 0, 0, 1 };
        
    mat_mul(mat, t, 4, mat);
    
}

void mat_rotx(float mat[16], float theta){
    mat_rotx2(mat, sin(theta), cos(theta));
}

void mat_rotx2(float mat[16], float sintheta, float costheta){
    float rx[16] = {
        1,        0,         0, 0,
        0, costheta, -sintheta, 0,
        0, sintheta,  costheta, 0,
        0,        0,         0, 1 };
    mat_mul(mat, rx, 4, mat);
}

void mat_roty(float mat[16], float theta){
    mat_roty2(mat, sin(theta), cos(theta));
}

void mat_roty2(float mat[16], float sintheta, float costheta){
    float ry[16] = {
       costheta, 0, sintheta, 0,
              0, 1,        0, 0,
      -sintheta, 0, costheta, 0,
              0, 0,        0, 1 };
    mat_mul(mat, ry, 4, mat);
}

void mat_rotz(float mat[16], float theta){
    mat_rotz2(mat, sin(theta), cos(theta));
}

void mat_rotz2(float mat[16], float sintheta, float costheta){
    float rz[16] = {
      costheta, -sintheta, 0, 0,
      sintheta,  costheta, 0, 0,
             0,         0, 1, 0,
             0,         0, 0, 1 };
    mat_mul(mat, rz, 4, mat);
    
}

void mat_scale(float mat[16], float sx, float sy, float sz){
    float s[16] = {
        sx,  0,  0, 0,
         0, sy,  0, 0,
         0,  0, sz, 0,
         0,  0,  0, 1 };
    mat_mul(mat, s, 4, mat);
    
}

void calc_viewport(float *xmin, float *ymin, float *xmax, float *ymax){
	for (int i = 0; i < N; ++i)
	{
		if(triangles[i].p1[0] < *xmin) *xmin = triangles[i].p1[0];
		if(triangles[i].p1[1] < *ymin) *ymin = triangles[i].p1[1];
		if(triangles[i].p1[0] > *xmax) *xmax = triangles[i].p1[0];
		if(triangles[i].p1[1] > *ymax) *ymax = triangles[i].p1[1];

		if(triangles[i].p2[0] < *xmin) *xmin = triangles[i].p2[0];
		if(triangles[i].p2[1] < *ymin) *ymin = triangles[i].p2[1];
		if(triangles[i].p2[0] > *xmax) *xmax = triangles[i].p2[0];
		if(triangles[i].p2[1] > *ymax) *ymax = triangles[i].p2[1];

		if(triangles[i].p3[0] < *xmin) *xmin = triangles[i].p3[0];
		if(triangles[i].p3[1] < *ymin) *ymin = triangles[i].p3[1];
		if(triangles[i].p3[0] > *xmax) *xmax = triangles[i].p3[0];
		if(triangles[i].p3[1] > *ymax) *ymax = triangles[i].p3[1];
	}
}

void calc_mat(float m[16]){
	mat_id(m);
	m[10] = 0; //proyectar en xy
	
	float xmin, ymin, xmax, ymax;
	calc_viewport(&xmin, &ymin, &xmax, &ymax);
	float w = (xmax-xmin);
	float h = (ymax-ymin);
	float sx = width/w, sy = height/h, dx = 0, dy = 0;
	if(sx < sy) sy = sx;
	else sx = sy;
	
	/*if(h > w){
		sx = sy*w/h;
		dx = width*(1-h/w)/2;
	} else {
		sy = sx*h/w;
		dy = height*(1-w/h)/2;
	}
	*/
	
	float dx = (width-(xmax/sx))/2, dy = ((ymax/sy)-height)/2;
	mat_trans(m, dx, dy, 0);
	mat_scale(m, sx, sy, 1);
	mat_trans(m, -xmin, -ymin, 0);
	for(int i=0; i<4; ++i)
		printf("[%f,%f,%f,%f]\n", m[4*i], m[4*i+1], m[4*i+2], m[4*i+3]);
	
}

void proyectar(float m[16], float a[], struct Trianglew t){
	float r1[4], r2[4], r3[4];
	mat_mul(m, t.p1, 1, r1);
	mat_mul(m, t.p2, 1, r2);
	mat_mul(m, t.p3, 1, r3);
	draw_pixel(a, r1[0], r1[1], 1,1,1);
	draw_pixel(a, r2[0], r2[1], 1,1,1);
	draw_pixel(a, r3[0], r3[1], 1,1,1);
	//dibujar_linea(a, (int)r1[0], (int)r1[1], (int)r2[0], (int)r2[1]);
	//dibujar_linea(a, (int)r1[0], (int)r1[1], (int)r3[0], (int)r3[1]);
	//dibujar_linea(a, (int)r3[0], (int)r3[1], (int)r2[0], (int)r2[1]);
}

void dibujar_linea(float a[], int _x1, int _y1, int _x2, int _y2){
	if( _x1 < 0 || _x1 > width || _x2 < 0 || _x2 > width ||
		_y1 < 0 || _y1 > height|| _y2 < 0 || _y2 > height )
		printf("(%d,%d), (%d,%d)\n", _x1, _y1, _x2, _y2);
	if(_x2-_x1 == 0){
		int ymin = _y1<_y2 ? _y1 : _y2;
		int ymax = _y1<_y2 ? _y2 : _y1;
		for(int i=ymin; i<=ymax; ++i){
			draw_pixel(a, _x1, i, 1, 1, 1);
		}
		//glDrawPixels(width, height, GL_RGB, GL_FLOAT, a);
		//glutSwapBuffers();
		return;
	}
	float m = (float)(_y2-_y1)/(_x2-_x1);
	if(m < -1)
		pendiente_menor_menos1(a, _x1, _y1, _x2, _y2);
	else if (m < 0)
		pendiente_entre_menos1y0(a, _x1, _y1, _x2, _y2);
	else if (m < 1)
		pendiente_entre_1y0(a, _x1, _y1, _x2, _y2);
	else pendiente_mayor_1(a, _x1, _y1, _x2, _y2);
}

void pendiente_menor_menos1(float *a, int x1, int y1, int x2, int y2){
	int dx = x2-x1, dy = y2-y1;
	int bdx = y1*dx - x1*dy;
	int x = x1<x2 ? x1 : x2;
	int y = x1<x2 ? y1 : y2;
	draw_pixel(a, x, y, 1, 1, 1);
	int ymax = x1<x2 ? y2 : y1;
	while(y>ymax){
		if(dy*(x+0.5f)-dx*(y-1)+bdx > 0)
			draw_pixel(a, ++x, --y, 1, 1, 1);
		else draw_pixel(a, x, --y, 1, 1, 1);
		
	}
}

void pendiente_entre_menos1y0(float *a, int x1, int y1, int x2, int y2){
	int dx = x2-x1, dy = y2-y1;
	int bdx = y1*dx - x1*dy;
	int x = x1<x2 ? x1 : x2;
	int y = x1<x2 ? y1 : y2;
	draw_pixel(a, x, y, 1, 1, 1);
	int xmax = x1<x2 ? x2 : x1;
	while(x<xmax){
		if(dy*(x+1)-dx*(y-0.5f)+bdx < 0)
			draw_pixel(a, ++x, --y, 1, 1, 1);
		else draw_pixel(a, ++x, y, 1, 1, 1);
	}
}

void pendiente_entre_1y0(float *a, int x1, int y1, int x2, int y2){
	int dx = x2-x1, dy = y2-y1;
	int bdx = y1*dx - x1*dy;
	int x = x1<x2 ? x1 : x2;
	int y = x1<x2 ? y1 : y2;
	draw_pixel(a, x, y, 1, 1, 1);
	int xmax = x1<x2 ? x2 : x1;
	while(x<xmax){
		if(dy*(x+1)-dx*(y+0.5f)+bdx < 0)
			draw_pixel(a, ++x, y, 1, 1, 1);
		else draw_pixel(a, ++x, ++y, 1, 1, 1);
		
	}
}

void pendiente_mayor_1(float *a, int x1, int y1, int x2, int y2){
	int dx = x2-x1, dy = y2-y1;
	int bdx = y1*dx - x1*dy;
	int x = x1<x2 ? x1 : x2;
	int y = x1<x2 ? y1 : y2;
	draw_pixel(a, x, y, 1, 1, 1);
	int ymax = x1<x2 ? y2 : y1;
	//int *i;
	//printf("%d", sizeof(float)*height*width*3);
	//scanf("%d", i);
	while(y<ymax){
		if(dy*(x+0.5f)-dx*(y+1)+bdx < 0)
			draw_pixel(a, ++x, ++y, 1, 1, 1);
		else draw_pixel(a, x, ++y, 1, 1, 1);
	}
}

void setup(){
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
}

void draw_pixel(float *a, int x, int y, float r, float g, float b){
	printf("(%d, %d)\n", x, y);
	int p = y*width*3+x*3;
	a[p] = r;
	a[p+1] = g;
	a[p+2] = b;
}

void abrir_archivo_y_procesa(struct Trianglew tris[]){
	FILE *file;
	float buff;
	file = fopen("OBJETOS-3D/itokawa_f0049152.tri", "r");
	int i=0;
	while(fscanf(file, "%f", &buff) != EOF){
		tris[i].p1[0] = buff;
		fscanf(file, "%f", &buff);
		tris[i].p1[1] = buff;
		fscanf(file, "%f", &buff);
		tris[i].p1[2] = buff;
		tris[i].p1[3] = 1;
		fscanf(file, "%f", &buff);
		tris[i].p2[0] = buff;
		fscanf(file, "%f", &buff);
		tris[i].p2[1] = buff;
		fscanf(file, "%f", &buff);
		tris[i].p2[2] = buff;
		tris[i].p2[3] = 1;
		fscanf(file, "%f", &buff);
		tris[i].p3[0] = buff;
		fscanf(file, "%f", &buff);
		tris[i].p3[1] = buff;
		fscanf(file, "%f", &buff);
		tris[i].p3[2] = buff;
		tris[i].p3[3] = 1;
		fscanf(file, "%f", &buff);
		++i;
	}
	fclose(file);
}
