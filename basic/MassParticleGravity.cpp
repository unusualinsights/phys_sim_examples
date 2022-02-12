#include <GL/freeglut.h>

#include <iostream>

// A data type representing a point (x, y)
struct Point2D {
  GLdouble x;
  GLdouble y;
};

// A data type representing a vector (x, y)
typedef Point2D Vector2D;

// Acceleration due to gravity on the Moon
// https://en.wikipedia.org/wiki/Gravitation_of_the_Moon
const Vector2D kGMoon{.x = 0.0, .y = -1.625};

// Returns how many seconds have elapsed since glutInit was called.
double GetCurrentTimeInSeconds() { return glutGet(GLUT_ELAPSED_TIME) / 1000.0; }

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

// Returns the x and y coordinates of the velocity of a Particle at the end of a
// time step of size |dt| if its current (before the time step) velocity is
// |current_velocity| and it has the given |mass| and experiences the specified
// |net_force|.
Vector2D MakeUpdatedVelocity(const Vector2D& current_velocity, double mass,
                             const Vector2D& net_force, double dt) {
  // Update velocity: symplectic Euler integration step 1
  return Vector2D{.x = current_velocity.x + dt * net_force.x / mass,
                  .y = current_velocity.y + dt * net_force.y / mass};
}

// Returns the (x, y) position of a particle at the current time + |dt| given
// its |current_position| at the current time and its |velocity| returned by
// MakeUpdatedVelocity.
Point2D MakeUpdatedPosition(const Point2D& current_position, double dt,
                            const Vector2D& velocity) {
  // Update position: symplectic Euler integration step 2
  //
  // Compute the Particle's new position at the new time, |t| + |dt|.
  // x_(t + dt) = x_(t) + dt * u(x_, y_, 0, t)
  // y_(t + dt) = y_(t) + dt * v(x_, y_, 0, t)
  return Point2D{.x = current_position.x + dt * velocity.x,
                 .y = current_position.y + dt * velocity.y};
}

// Draws this Particle at (current_position.x, current_position.y, 0.0).
void DrawParticleAt(const Point2D& current_position) {
  SetColorToRed();
  DrawSphereCenteredAt(current_position.x, current_position.y, 0.0);
}

// A data type representing a particle (x(t), y(t), 0)
//
// An instance of this class, i.e., a single Particle object in the computer's
// memory, will draw itself onto the screen whenever its Render() method is
// called.
class Particle {
 public:
  // Creates a Particle with the given |mass| starting at position |pos0| with
  // initial velocity |vel0|.
  Particle(double mass, Point2D pos0, Vector2D vel0)
      : mass_(mass),
        weight_{.x = kGMoon.x * mass, .y = kGMoon.y * mass},
        pos_(pos0),
        vel_(vel0) {}

  // Computes this Particle's updated position and renders it onto the screen.
  void UpdateAndRender(double dt) {
    // Use the force acting on the Particle to update its velocity.
    //
    // Pass the net force acting on this Particle as the |net_force| parameter.
    // In this case, the only force acting on this Particle is its |weight_|.
    vel_ = MakeUpdatedVelocity(vel_, mass_, weight_, dt);

    // Use the updated velocity to update the Particle's position.
    pos_ = MakeUpdatedPosition(pos_, dt, vel_);

    DrawParticleAt(pos_);
  }

 private:
  const double mass_;      // This Particle's mass
  const Vector2D weight_;  // This Particle's weight on the Moon
  Point2D pos_;   // This Particle's current position, (pos_.x, pos_.y, 0.0)
  Vector2D vel_;  // This Particle's current velocity, (vel_.x, vel_.y, 0.0)
};

// A data type encapsulating (holding, owning) a Particle to simulate and a
// current time in the simulation
//
// This class carries a simulation forward in time and directs any Particles it
// owns to calculate their own positions at new time values and draw themselves
// at those positions as time elapses.
class ParticleSimulator {
 public:
  // Creates a Particle and a time value.
  ParticleSimulator()
      : particle_(0.01, Point2D{.x = -0.5, .y = -0.5},
                  Vector2D{.x = 0.5, .y = 2.0}),
        current_time_(GetCurrentTimeInSeconds()) {}

  // Updates the simulation time and directs any Particles to update and render
  // themselves at the new simulation time.
  void UpdateAndRender() {
    // Update dt to be the current time minus the old time.
    double next_time = GetCurrentTimeInSeconds();  // t + dt
    double dt = next_time - current_time_;         // dt = (t + dt) - t

    // Compute the Particle's position at |next_time| and draw it at that
    // location.
    particle_.UpdateAndRender(dt);

    // It's the new time now, but it will be the old time next time this method
    // is called!
    current_time_ = next_time;
  }

 private:
  Particle particle_;    // The particle whose motion we're simulating
  double current_time_;  // The last time the particle's position was updated
};

// Global pointer to the object running the particle simulation
ParticleSimulator* particle_simulator;

// Tells OpenGL exactly what to redraw in each frame.
void RenderScene() {
  if (!particle_simulator) {
    std::cout << "ERROR: particle_simulator is not initialized!" << std::endl;
    return;
  }

  DrawBackground();
  TurnOnLight();

  particle_simulator->UpdateAndRender();

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

  int win = glutCreateWindow("Particle Flying over the Moon");
  std::cout << "Window with ID " << win << " opened successfully." << std::endl;

  particle_simulator = new ParticleSimulator();

  glutDisplayFunc(RenderScene);
  glutIdleFunc(RenderScene);

  glutMainLoop();

  return 0;
}
