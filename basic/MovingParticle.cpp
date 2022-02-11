#include <GL/freeglut.h>

#include <cmath>
#include <iostream>

// Sets the background color to black.
void SetBackgroundColor() {
  const GLclampf kClearRed = 0.0f;
  const GLclampf kClearGreen = 0.0f;
  const GLclampf kClearBlue = 0.0f;
  const GLclampf kClearAlpha = 1.0f;
  glClearColor(kClearRed, kClearGreen, kClearBlue, kClearAlpha);
}

// Sets the background color and a few other display settings for each new frame
// OpenGL draws.
void DrawBackground() {
  SetBackgroundColor();

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
}

// Enables one light source.
void TurnOnLight() {
  const GLfloat kAmbient[4] = {0.2f, 0.2f, 0.2f, 1.0f};
  const GLfloat kDiffuse[4] = {0.4f, 0.4f, 0.4f, 1.0f};
  const GLfloat kSpecular[4] = {0.7f, 0.7f, 0.7f, 1.0f};
  const GLfloat kPosition[4] = {1.0f, 1.0f, 0.5f, 1.0f};
  glLightfv(GL_LIGHT0, GL_AMBIENT, kAmbient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, kDiffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, kSpecular);
  glLightfv(GL_LIGHT0, GL_POSITION, kPosition);
  glEnable(GL_LIGHT0);
}

// Sets the ambient, diffuse, and specular RGBA reflectances of the material
// model applied to any vertices drawn subsequently.
void SetMaterialColor(const GLfloat material_color[4]) {
  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, material_color);
  glMaterialfv(GL_FRONT, GL_SPECULAR, material_color);
}

// Sets the color of any vertices we draw to red.
void SetColorToRed() {
  const GLfloat kRedColor[4] = {1.0f, 0.0f, 0.0f, 1.0f};
  SetMaterialColor(kRedColor);
}

// Draws a sphere whose center is at (center_x, center_y, center_z) in the 3D
// coordinate system in which we're drawing geometry.
void DrawSphereCenteredAt(GLdouble center_x, GLdouble center_y,
                          GLdouble center_z) {
  glPushMatrix();
  glTranslated(center_x, center_y, center_z);

  const GLdouble kRadius = 0.1;
  const GLint kLongitudeSlices = 50;
  const GLint kLatitudeStacks = 50;
  glutSolidSphere(kRadius, kLongitudeSlices, kLatitudeStacks);

  glPopMatrix();
}

// Tells OpenGL exactly what to redraw in each frame.
void RenderScene() {
  DrawBackground();
  TurnOnLight();

  // t = # of seconds since glutInit was called
  double t = glutGet(GLUT_ELAPSED_TIME) / 1000.0;

  // Compute a new position for the sphere that will move it around in a circle
  // of radius 0.5 about the origin as time, |t|, elapses.
  GLdouble x = 0.5 * cos(t);
  GLdouble y = 0.5 * sin(t);

  SetColorToRed();

  // Draw the sphere centered at this new position, (x, y).
  DrawSphereCenteredAt(x, y, 0.0);

  glutSwapBuffers();
}

int main(int argc, char** argv) {
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);

  const int kWidth = 500;
  const int kHeight = 500;
  glutInitWindowSize(kWidth, kHeight);

  const int kTopLeftX = 200;
  const int kTopLeftY = 100;
  glutInitWindowPosition(kTopLeftX, kTopLeftY);

  int win = glutCreateWindow("Moving Particle");
  std::cout << "Window with ID " << win << " opened successfully." << std::endl;

  glutDisplayFunc(RenderScene);
  glutIdleFunc(RenderScene);

  glutMainLoop();

  return 0;
}
