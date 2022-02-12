#include <GL/freeglut.h>

#include <iostream>
#include <vector>

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

// Draws a sphere with radius |radius| whose center is at (center_x, center_y,
// center_z).
void DrawSphereCenteredAt(GLdouble center_x, GLdouble center_y,
                          GLdouble center_z, GLdouble radius) {
  glPushMatrix();
  glTranslated(center_x, center_y, center_z);

  const GLint kLongitudeSlices = 50;
  const GLint kLatitudeStacks = 50;
  glutSolidSphere(radius, kLongitudeSlices, kLatitudeStacks);

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

// Returns the (x, y) position of a particle at time |current_time| + |dt| given
// its |current_position| at |current_time| and its |velocity| returned by
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

// Draws this Particle at (current_position.x, current_position.y, 0.0) with
// radius |radius|.
void DrawParticleAt(const Point2D& current_position, double radius) {
  SetColorToRed();
  DrawSphereCenteredAt(current_position.x, current_position.y, 0.0, radius);
}

// A data type representing a particle (x(t), y(t), 0)
//
// An instance of this class, i.e., a single Particle object in the computer's
// memory, will draw itself onto the screen whenever its Render() method is
// called.
class Particle {
 public:
  // Creates a Particle with the given |mass| starting at position |pos0| with
  // initial velocity |vel0|. When displayed, this Particle will have the
  // specified |radius|.
  Particle(double mass, double radius, Point2D pos0, Vector2D vel0)
      : mass_(mass),
        weight_{.x = kGMoon.x * mass, .y = kGMoon.y * mass},
        radius_(radius),
        pos_(pos0),
        vel_(vel0) {}

  double mass() const { return mass_; }
  double radius() const { return radius_; }
  double x() const { return pos_.x; }
  double y() const { return pos_.y; }
  double u() const { return vel_.x; }
  double v() const { return vel_.y; }

  // Clears this Particle's force accumulator.
  void ZeroOutForceAccumulator() { force_accumulator_ = {.x = 0.0, .y = 0.0}; }

  // Adds the force |f| to this Particle's force accumulator.
  void ApplyForce(const Vector2D& f) {
    force_accumulator_.x += f.x;
    force_accumulator_.y += f.y;
  }

  // Computes this Particle's updated position and renders it onto the screen.
  //
  // Make sure this Particle's force accumulator is zeroed out and updated with
  // the forces acting on it in the current time step before calling this
  // method.
  void UpdatePosition(double dt) {
    // Use the net force acting on the Particle to update its velocity.
    //
    // So, pass the |force_accumulator_| as the |net_force| parameter to this
    // function.
    vel_ = MakeUpdatedVelocity(vel_, mass_, force_accumulator_, dt);

    // Use the updated velocity to update the Particle's position.
    pos_ = MakeUpdatedPosition(pos_, dt, vel_);
  }

  // Renders this Particle.
  void Render() { DrawParticleAt(pos_, radius_); }

 private:
  const double mass_;      // This Particle's mass
  const Vector2D weight_;  // This Particle's weight on the Moon
  const GLdouble radius_;  // This Particle's drawn sphere's radius
  Point2D pos_;   // This Particle's current position, (pos_.x, pos_.y, 0.0)
  Vector2D vel_;  // This Particle's current velocity, (vel_.x, vel_.y, 0.0)

  // Sum of forces acting on this Particle
  //
  // This force accumulator will be cleared at the beginning of each time step.
  // Then as forces acting on this Particle are calculated, they will be added
  // to this 2D vector to be used in updating this Particle's velocity and
  // position at the end of that time step.
  Vector2D force_accumulator_;
};

// Prints info of |p|, the |particle_number|th Particle created for simulation.
void Print(std::size_t particle_number, const Particle& p) {
  std::cout << std::endl;
  std::cout << "Particle #" << particle_number << ":" << std::endl;
  std::cout << "- (x0, y0) = (" << p.x() << ", " << p.y() << ")" << std::endl;
  std::cout << "- (u0, v0) = (" << p.u() << ", " << p.v() << ")" << std::endl;
  std::cout << "- mass = " << p.mass() << std::endl;
  std::cout << "- radius = " << p.radius() << std::endl;
  std::cout << std::endl;
}

// Initializes and returns Particles to begin the particle simulation.
std::vector<Particle> MakeInitialParticles() {
  const std::size_t kNumParticles = 3;

  std::vector<Particle> particles;

  if (kNumParticles < 2) {
    std::cout << "ERROR: Need >= 2 Particles!" << std::endl;
    return particles;  // return empty list of Particles
  }

  const double kInitPosSpacing = 2.0 / (kNumParticles + 1);
  const double kSpeedSpacing = 1.0 / (kNumParticles - 1);
  for (std::size_t p_index = 0; p_index < kNumParticles; p_index++) {
    double x0 = kInitPosSpacing * (p_index + 1) - 1.0;
    Point2D pos0{.x = x0, .y = -0.5};

    double speed_interp_ratio = p_index * kSpeedSpacing;

    double v0 = 1.5 + speed_interp_ratio * 0.5;
    Vector2D vel0{.x = -x0, .y = v0};

    double mass = 0.01 + speed_interp_ratio * 0.01;
    double radius = 0.1 + speed_interp_ratio * 0.05;

    Particle p(mass, radius, pos0, vel0);
    particles.push_back(p);

    Print(p_index + 1, p);
  }

  std::cout << particles.size() << " particles created." << std::endl;
  std::cout << std::endl;

  return particles;
}

// Returns whether
// ||(q.x(), q.y()) - (p.x(), p.y())|| < p.radius() + q.radius(),
// i.e., whether Particles |p| and |q| are colliding.
bool Colliding(const Particle& p, const Particle& q) {
  double dx = q.x() - p.x();
  double dy = q.y() - p.y();
  double dist_squared_btw_centers = dx * dx + dy * dy;
  double sum_of_radii = p.radius() + q.radius();
  return dist_squared_btw_centers < sum_of_radii * sum_of_radii;
}

// Returns the average impulse force resulting from the change in momentum from
// mass * old velocity to mass * new velocity, i.e., from
// me.mass() * (me.u(), me.v()) to me.mass() * (new_u, new_v), assuming the
// force is applied during a duration |dt|.
Vector2D MakeAvgCollisionForce(double new_u, double new_v, const Particle& me,
                               double dt) {
  // Impulse-Momentum Theorem: Impulse equals change in momentum, which equals
  // change in velocity times mass
  double impulse_x = (new_u - me.u()) * me.mass();
  double impulse_y = (new_v - me.v()) * me.mass();

  // Impulse = average collision force * dt, so average collision force during
  // the time increment, dt, is impulse / dt.
  return Vector2D{.x = impulse_x / dt, .y = impulse_y / dt};
}

// Returns the action, i.e., the force acting on |me|, resulting from an elastic
// collision with |other|.
//
// We use the formulas for velocities after impact in terms of the particle
// masses and pre-impact velocities and positions from the section,
// "Two-dimensional collision with two moving objects" from
// https://en.wikipedia.org/wiki/Elastic_collision.
Vector2D MakeCollisionForce(const Particle& me, const Particle& other,
                            double dt) {
  double mass_factor = 2.0 * other.mass() / (me.mass() + other.mass());
  double du = me.u() - other.u();
  double dv = me.v() - other.v();
  double dx = me.x() - other.x();
  double dy = me.y() - other.y();
  double d_vel_dot_d_pos = du * dx + dv * dy;
  double vel_pos_factor = d_vel_dot_d_pos / (dx * dx + dy * dy);
  double factor = mass_factor * vel_pos_factor;
  double new_u = me.u() - factor * dx;
  double new_v = me.v() - factor * dy;
  return MakeAvgCollisionForce(new_u, new_v, me, dt);
}

// Returns whether the Particle |me| hit a side wall of the imaginary tunnel
// within which we will force the particles to stay in this simulation.
bool HitSideWall(const Particle& me) {
  const double kLeftWallX = -1.0;
  const double kRightWallX = 1.0;

  return me.x() - me.radius() < kLeftWallX ||
         me.x() + me.radius() > kRightWallX;
}

// Returns whether the Particle |me| hit the top or bottom wall of the imaginary
// tunnel within which we will force the particles to stay in this simulation.
bool HitTopOrBottomWall(const Particle& me) {
  const double kBottomWallY = -1.0;
  const double kTopWallY = 1.0;

  return me.y() - me.radius() < kBottomWallY ||
         me.y() + me.radius() > kTopWallY;
}

// Returns whether the Particle |me| has hit any wall of the imaginary tunnel
// within which we will force the particles to stay in this simulation.
bool HitWall(const Particle& me) {
  return HitSideWall(me) || HitTopOrBottomWall(me);
}

// Returns the force applied by a wall of the viewing frustum (the imaginary
// tunnel represented by the window displaying our simulation) onto the Particle
// |me| due to collision with the wall.
Vector2D MakeWallCollisionForce(const Particle& me, double dt) {
  double new_u = me.u();
  // Note we assume the radius of the ball is less than both the width and
  // height of the viewing frustum.
  if (HitSideWall(me)) {
    new_u = -new_u;
  }

  double new_v = me.v();
  if (HitTopOrBottomWall(me)) {
    new_v = -new_v;
  }

  return MakeAvgCollisionForce(new_u, new_v, me, dt);
}

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
      : particles_(MakeInitialParticles()),
        current_time_(GetCurrentTimeInSeconds()) {}

  // Updates the simulation time and directs any Particles to update and render
  // themselves at the new simulation time.
  void UpdateAndRender() {
    // Update dt to be the current time minus the old time.
    double next_time = GetCurrentTimeInSeconds();  // t + dt
    double dt = next_time - current_time_;         // dt = (t + dt) - t

    UpdateParticleForces(dt);
    UpdateParticlePositions(dt);

    for (Particle& particle : particles_) {
      particle.Render();
    }

    // It's the new time now, but it will be the old time next time this method
    // is called!
    current_time_ = next_time;
  }

 private:
  // Updates the forces acting on each Particle in the simulation for a given
  // time step.
  void UpdateParticleForces(double dt) {
    // Clear all forces acting on all Particles.
    for (Particle& particle : particles_) {
      particle.ZeroOutForceAccumulator();
    }

    // Apply Particle-to-Particle collision forces.
    for (std::size_t p_index = 0; p_index < particles_.size(); p_index++) {
      Particle& particle = particles_.at(p_index);

      // If |particle| is colliding with any other Particle, then apply the
      // appropriate action force to |particle| and reaction force (negation of
      // the action force) to the other Particle.
      //
      // Start the loop searching for other Particles at p_index + 1 so we visit
      // a each pair of Particles exactly once.
      for (std::size_t q_index = p_index + 1; q_index < particles_.size();
           q_index++) {
        Particle& other_particle = particles_.at(q_index);

        if (!Colliding(particle, other_particle)) {
          continue;
        }

        // |particle| is colliding with |other_particle|.
        // So, apply the action force to |particle| ...
        Vector2D action_force =
            MakeCollisionForce(particle, other_particle, dt);
        particle.ApplyForce(action_force);

        // ... and the reaction force to |other_particle|.
        Vector2D reaction_force{.x = -action_force.x, .y = -action_force.y};
        other_particle.ApplyForce(reaction_force);
      }
    }

    // Apply Particle-wall collision forces.
    for (Particle& particle : particles_) {
      if (!HitWall(particle)) {
        continue;
      }

      // |particle| hit the wall, so apply a wall-to-Particle collision force.
      //
      // We don't apply a reaction force from |particle| to the wall since we
      // assume the wall has infinite mass so any finite force applied to it
      // results in an acceleration of zero for the wall itself.
      Vector2D wall_force = MakeWallCollisionForce(particle, dt);
      particle.ApplyForce(wall_force);
    }
  }

  // Directs each Particle in the simulation to update its position based on the
  // forces it is experiencing.
  void UpdateParticlePositions(double dt) {
    for (Particle& particle : particles_) {
      particle.UpdatePosition(dt);
    }
  }

  std::vector<Particle> particles_;  // The Particles we're simulating
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
