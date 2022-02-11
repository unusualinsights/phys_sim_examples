#include <GL/freeglut.h>

#include <iostream>

// Sets the background color to black.
//
// OpenGL keeps redrawing our graphics in a non-stop loop. Before redrawing each
// frame of our scene, we need to tell OpenGL what color we want the background
// to be.
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

  // Tell OpenGL exactly what to clear before drawing each frame.
  //
  // |GL_COLOR_BUFFER_BIT| tells OpenGL to clear the colors of pixels in tts
  // display window to the color set by the above glClearColor function call.
  //
  // |GL_DEPTH_BUFFER_BIT| tells OpenGL to clear any depth information about the
  // geometry it's currently displaying. This makes sure OpenGL doesn't block
  // your view of geometry it draws in future frames with geometry that was
  // closest to you in previous frames.
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Tell OpenGL to draw geometry that's closer to you, the viewer, in front of
  // geometry that's farther away from you.
  glEnable(GL_DEPTH_TEST);

  // Tell OpenGL not to waste time drawing geometry we can't see if it's blocked
  // by other geometry that's closer to us.
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
}

// Sets the color of any geometry we draw to red.
void SetColorToRed() {
  const GLfloat kRedColor[4] = {1.0f, 0.0f, 0.0f, 1.0f};
  glColor4fv(kRedColor);
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

  SetColorToRed();
  DrawSphereCenteredAtOrigin();

  // In the call to glutInitDisplayMode below, we use GLUT to tell OpenGL to use
  // "double buffering." This means OpenGL will keep two complete copies of
  // pixels: a "front buffer" of pixels currently displayed on the screen and a
  // "back buffer" where it sets pixel values you can't see on the screen. When
  // glutSwapBuffers() is called, the back and front buffers switch places,
  // causing whatever was in the old back buffer (which is now the front buffer)
  // to be displayed on the screen.
  glutSwapBuffers();
}

int main(int argc, char** argv) {
  // Make sure GLUT is ready to work with the computer's windowing system!
  glutInit(&argc, argv);

  // Set display mode settings.
  //
  // |GLUT_DEPTH| tells GLUT we want closer objects hiding farther away objects.
  //
  // |GLUT_DOUBLE| tells GLUT to use double buffering. See the above call to
  // glutSwapBuffers().
  //
  // |GLUT_RGBA| tells GLUT to use red, green, and blue color channels plus an
  // "alpha" channel to allow geometry to be opaque, translucent, or
  // transparent.
  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);

  // Set the size of the window this program will open, in pixels.
  const int kWidth = 500;
  const int kHeight = 500;
  glutInitWindowSize(kWidth, kHeight);

  // Set the coordinates, in pixels, of the top left corner of the window this
  // program will open. Note that (0, 0) is at the top left corner of your
  // computer screen and the y-coordinate increases *downward* on your screen.
  // So, we'll end up with a window whose top left corner is 200 pixels to the
  // right of and 100 pixels below the top left corner of the computer screen.
  const int kTopLeftX = 200;
  const int kTopLeftY = 100;
  glutInitWindowPosition(kTopLeftX, kTopLeftY);

  // Create the window with the specified title.
  int win = glutCreateWindow("Hello, Sphere!");
  std::cout << "Window with ID " << win << " opened successfully." << std::endl;

  // Tell GLUT what method to call to redraw the scene in each frame.
  glutDisplayFunc(RenderScene);
  glutIdleFunc(RenderScene);

  // Keep redrawing the scene until this program is terminated.
  glutMainLoop();

  // Everything's a-okay!
  return 0;
}
