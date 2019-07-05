#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <math.h> 


using std::ostream;
using std::istream;
using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::string;

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
class vector2d {
  public:
  double x,y;

  vector2d() {
   x=y=0;
  }
  vector2d(double xnum, double ynum) {
    x = xnum;
    y = ynum;
  }

  vector2d operator-(const vector2d& other) const {
    return vector2d(x-other.x,y-other.y);
  }
 friend vector2d operator*(const double &num, const vector2d &v) { //looks good
    return vector2d(num*v.x, num*v.y);
  }

  friend double operator*(const vector2d &u, const vector2d &v) { //confirmed
    return u.x*v.x + u.y*v.y;
  }
};

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
class Stone {
  public:
    vector2d pos, vel;
    double mass;
    double radius;
    string name;

    vector2d momentum() const {
        vector2d momentum;
        momentum.x = vel.x * mass;
        momentum.y = vel.y * mass;
        return momentum;
    }

    double energy() const {
        return (vel.x * vel.x + vel.y * vel.y) * mass * 0.5;
    }

    // create a stone according to input from cin
    void createstone(double ma, double rad, double px, double py,
    double vx, double vy, string na) {
        mass = ma;
        radius = rad;
        pos.x = px;
        pos.y = py;
        vel.x = vx;
        vel.y = vy;
        name = na;
        return;
    }

    void show() {
        cout << name << " m=" << mass << " R=" << radius << " p=(" <<
        pos.x << "," << pos.y << ") v=(" << vel.x << "," << vel.y <<
        ")" << endl;
    }

    // move the stone to the collisiion time point
    void move(double time) {
        pos.x += vel.x*time;
        pos.y += vel.y*time;
    }
};

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
class Collision {
  public:
    Stone one;
    Stone two;
    double time;

    bool valid() {
        // calculates the position between the two stones
        vector2d u;
        u.x = one.pos.x - two.pos.x;
        u.y = one.pos.y - two.pos.y;
        //calculates the velocity between the two stones
        vector2d v;
        v.x = one.vel.x - two.vel.x;
        v.y = one.vel.y - two.vel.y;

        double check = u.x*v.x + u.y*v.y;
        if (check < 0.00) {
            return true;
        }
        return false;
    } 

    // define
    Collision () {
        time = 0;
    }

    Collision(double t, Stone s1, Stone s2) {
        time = t;
        one = s1;
        two = s2;
  }
};
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// a class for the vector of stones and updates after collision
class StoneList {
  public:
    vector<Stone> stones;
    StoneList(vector<Stone> ston) {
        stones  = ston;
    }

    void collide(Stone one, Stone two) {
        vector2d vi = one.vel;
        vector2d vj = two.vel;

        double m1factor = (2*two.mass)/(one.mass+two.mass);
        double m2factor = (2*one.mass)/(one.mass+two.mass);

        vector2d w1;
        w1.x = one.vel.x - two.vel.x;
        w1.y = one.vel.y - two.vel.y;
        vector2d w2;
        w2.x = two.vel.x - one.vel.x;
        w2.y = two.vel.y - one.vel.y;

        vector2d u1;
        u1.x = one.pos.x - two.pos.x;
        u1.y = one.pos.y - two.pos.y;

        vector2d u2;
        u2.x = two.pos.x - one.pos.x;
        u2.y = two.pos.y - one.pos.y;

        // calculate new velocitiies
        vector2d via = one.vel - ((m1factor)*((w1*u1)/(u1*u1))*u1);

        vector2d vja = two.vel - ((m2factor)*((w2*u2)/(u2*u2))*u2);

        one.vel = via;
        two.vel = vja;

        // use string name to detect which stone to update
        vector<string> namestring;
        for (int x = 0; x < stones.size(); x++) {
            namestring.push_back(stones.at(x).name);
        }
        for (int x = 0; x < stones.size(); x++) {
            for (int y = 0; y < namestring.size(); y++) {
                if (namestring.at(x) == one.name) {
                    stones.at(x).vel = one.vel;
                }
                if (namestring.at(x) == two.name) {
                    stones.at(x).vel = two.vel;
                }
            }
        }
    }
};
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// calculate when the stones will crash
double collide_time(Stone one, Stone two) {
    vector2d u;
    u.x = one.pos.x-two.pos.x;
    u.y = one.pos.y-two.pos.y;
    vector2d w;
    w.x = one.vel.x - two.vel.x;
    w.y = one.vel.y - two.vel.y;
    double rad = one.radius + two.radius;
    double radsq = (rad)*(rad);
    double c = (u.x * u.x + u.y * u.y) - radsq;
    double b = 2*(w.x*u.x + w.y *u.y);
    double a = (w.x*w.x+w.y*w.y);
    double discr = (b*b) - (4*a*c);
    double pos  = (-b + sqrt(discr))/(2*a);
    double neg  = (-b - sqrt(discr))/(2*a);

    if (neg > 0.0) {
      return neg;
    }
    if (pos <= 0.0) { 
        //if there is no collision, break the loop in main loop
        return 1001;
    }
    return pos;
}
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
Collision get_next_collision(vector<Stone> stones) {
    double time;
    double first = 1000.0;

    Stone one;
    Stone two;

    for (int i = 0; i < stones.size(); i++) {
        for (int j = i+1; j < stones.size(); j++) {
            Stone s1 = stones.at(i);
            Stone s2 = stones.at(j);
            time = collide_time(s1,s2);
            if (time < first) {
                first = time;
                 one = s1;
                two = s2;
            }
            }
    }
    // create collision event
    Collision c(first, one, two);
    return c;
}
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void show_stones(vector<Stone> stones) {
  double energy = 0;
  vector2d momentum;
  for (int x = 0; x < stones.size(); x++) {
      stones.at(x).show();
      momentum.x = momentum.x + stones.at(x).momentum().x;
      momentum.y = momentum.y + stones.at(x).momentum().y;
      energy += stones.at(x).energy();
  }
  cout << "energy: " << energy << "\n";
  cout << "momentum: (" << momentum.x << "," << momentum.y << ")\n";
}

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// sort the vector of stones into alphabetical order
vector<Stone> sortvector(vector<Stone> stones) {
    vector<Stone> newstone;
    vector<string> forsort;
    for (int x = 0; x < stones.size(); x++) {
        forsort.push_back(stones.at(x).name);
    }
    sort(forsort.begin(),forsort.end());
    for (int x = 0; x < forsort.size(); x++) {
        for (int y = 0; y < stones.size(); y++) {
            if (forsort.at(x) == stones.at(y).name) {
                newstone.push_back(stones.at(y));
            }
        }
    }
    return newstone;
}

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int main() {
  double overall_time = 0;
  double mass, rad, posx, posy, velx, vely;
  string name;
  double momentum = 0, energy = 0;
  Stone st;

  cout << "Please enter the mass, radius, x/y position, x/y velocity\n";
  cout << "and name of each stone\n";
  cout << "When complete, use EOF to stop entering\n";
  cout << "Ctrl + D will not work\n";

  // save input stuff into a vector
  // The Ctrl + D will cause the vector of stones to push in extra times
  // Can type chars to stop input
  vector<Stone> stones;
  while (cin) {
    cin >> mass >> rad >> posx >> posy >> velx >> vely >> name;
    if (mass > 0) {
        st.createstone(mass, rad, posx, posy, velx, vely, name);
        stones.push_back(st);
    } else {
        break;
    }
  }
  stones = sortvector(stones);
  // convert into class for later use
  StoneList stonevec = stones;
  cout << "\nHere are the initial stones.\n" <<  endl;
  show_stones(stonevec.stones);
  cout << "\nHere are the collision events.\n";

  // loop that check collisiion and print info
  while (true) {
    Collision c = get_next_collision(stonevec.stones);
    if (!c.valid()) {
        break;
    }
    if (c.time >= 1000) {
        break;
    } 
    double newtime = c.time;
    for (auto & st : stonevec.stones) {
        st.move(newtime);
    }
    overall_time += newtime;
    cout << "\ntime of event: " << overall_time << "\n";
    cout << "colliding " << c.one.name << " " << c.two.name << endl;

    stonevec.collide(c.one,c.two);

    show_stones(stonevec.stones);
  }
  return 1;
}

/* Sample Run of file:

Please enter the mass, radius, x/y position, x/y velocity
and name of each stone
When complete, use EOF to stop entering
Ctrl + D will not work
10 1 0 0 1 0 one
1 1 20 0 0 0 two
30 2 10 0 0.2 0 three
EOF

Here are the initial stones.

one m=10 R=1 p=(0,0) v=(1,0)
three m=30 R=2 p=(10,0) v=(0.2,0)
two m=1 R=1 p=(20,0) v=(0,0)
energy: 5.6
momentum: (16,0)

Here are the collision events.

time of event: 8.75
colliding one three
one m=10 R=1 p=(8.75,0) v=(-0.2,0)
three m=30 R=2 p=(11.75,0) v=(0.6,0)
two m=1 R=1 p=(20,0) v=(0,0)
energy: 5.6
momentum: (16,0)

time of event: 17.5
colliding three two
one m=10 R=1 p=(7,0) v=(-0.2,0)
three m=30 R=2 p=(17,0) v=(0.56129,0)
two m=1 R=1 p=(20,0) v=(1.16129,0)
energy: 5.6
momentum: (16,0) */
