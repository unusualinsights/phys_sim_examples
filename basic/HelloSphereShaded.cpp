#include <GL/freeglut.h>

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
  glEnable(GL_LIGHTING);  // allow turning on lights for some 3D shading!
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
}

// Enables one light source.
//
// In RGBA mode, the lighted color of a vertex is:
//
// Color = Material emission intensity +
//         Material ambient reflectance *
//           Lighting model full-scene ambient +
//         Contribution of each enabled light source
//
// Contribution of light source = Ambient + Diffuse + Specular
// Ambient = Material ambient reflectance * Light's ambient intensity
// Diffuse = Material diffuse reflectance * Light's diffuse intensity *
//           DotProduct(Vertex normal, Normalized vector from vertex to light)
// Specular = Material specular reflectance * Light's specular intensity *
//            (DotProduct(Normalized vertex-to-eye vector,
//                        Normaliezd vertex-to-light vector))^material_shininess
// Each of these 3 light source components attenuate equally based on:
// - Distance from vertex to light source
// - Light source direction
// - Spread exponent
// - Spread cutoff angle
// All dot products are clamped to zero if they're negative.
// Alpha value = alpha value of material diffuse reflectance
//
// https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glLightModel.xml
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

// Draws a sphere whose center is at the origin of the 3D coordinate system in
// which we're drawing geometry.
void DrawSphereCenteredAtOrigin() {
  const GLdouble kRadius = 0.25;
  const GLint kLongitudeSlices = 50;
  const GLint kLatitudeStacks = 50;
  glutSolidSphere(kRadius, kLongitudeSlices, kLatitudeStacks);
}

// Tells OpenGL exactly what to redraw in each frame.
void RenderScene() {
  DrawBackground();
  TurnOnLight();

  SetColorToRed();
  DrawSphereCenteredAtOrigin();

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

  int win = glutCreateWindow("Hello, Sphere!");
  std::cout << "Window with ID " << win << " opened successfully." << std::endl;

  glutDisplayFunc(RenderScene);
  glutIdleFunc(RenderScene);

  glutMainLoop();

  return 0;
}
