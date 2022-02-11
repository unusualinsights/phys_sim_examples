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

  // Set the initial position of the sphere's center.
  // The 'static' keyword ensures these values are only assigned once the very
  // first time these lines of code are executed. But, the 'static' keyword
  // also ensures the values of these variables remain in memory between calls
  // to this method. So, |x| and |y| will always store the current position of
  // our particle (the center of the sphere we draw) at the current time, |t|.
  static GLdouble x = 0.0;
  static GLdouble y = 0.0;

  // Set initial t to be # of seconds since glutInit was called.
  // The 'static' keyword ensures |t| is only assigned this initial value once
  // the very first time this line of code is executed. But, the 'static'
  // keyword also ensures this value remains in memory between calls to this
  // method. So, effectively, this variable will always store the current time,
  // which we call |t|, of this simulation.
  static double t = glutGet(GLUT_ELAPSED_TIME) / 1000.0;

  // Update dt.
  // You're dealing with variable time steps (|dt| can change each time
  // RenderScene() is called) already! What a superstar you are!
  double tPlusDt = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
  double dt = tPlusDt - t;

  // Compute the new velocity vector (u(x, y, z, t), v(x, y, z, t)) at time |t|.
  // In case you're wondering, z is 0 here (note the z-coordinate we pass to the
  // method call below to draw the sphere).
  GLdouble u = -(x - 0.5 * cos(t * t) * cos(t));
  GLdouble v = -(y - 0.5 * cos(t * t) * sin(t));

  // Compute the new position of the sphere now that it has a new velocity!
  // This is the 'Forward Euler' integration scheme: your new position equals
  // your position plus how far you moved from the old to the new time, moving
  // at the speed you had at the old time.
  x += dt * u;  // x(t + dt) = x(t) + dt * u(x, y, z, t)
  y += dt * v;  // y(t + dt) = y(t) + dt * v(x, y, z, t)

  // Save this new time for the next time step.
  // This is |tPlusDt| now but next time this method is called, this will be the
  // new |t|.
  t = tPlusDt;

  SetColorToRed();
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

  int win = glutCreateWindow("Particle Moving in a Velocity Field");
  std::cout << "Window with ID " << win << " opened successfully." << std::endl;

  glutDisplayFunc(RenderScene);
  glutIdleFunc(RenderScene);

  glutMainLoop();

  return 0;
}
