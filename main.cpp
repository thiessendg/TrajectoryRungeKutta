#include <iostream>	// Standard Input Output
#include <cmath>	// Use Of Math Functions
#include <fstream>	// File Stream Input Output
#include <sstream>	// Used For Variable String Name
#include <iomanip>

using namespace std;

//constants
const double g = -9.80665;	//gravitational constant, the y acceleration
const double a = 0.0;		// x acceleration
const double pi = 3.14159265358979323846;
const double re = 6371000.0;//mean Radius of Earth in meters
const double deg2rad = pi / 180.0;

struct State	//the state of the projectile - position, velocity, and acceleration
{
	double xPos, yPos;
	double xVel, yVel;
	double xAcc, yAcc;
};

//function prototypes
void rk4(State & object, const double);
double gravity(const double);

//the main program
int main(int argc, char * argv[])
{
    //the input variables
    double initAlt;	//initial altitude
	double initVel;	//initial velocity at time0
	double theta;	//firing angle in degrees
	double timeStep;//the time step
	double duration;//time to run sim

	if (argc != 6) //if we got no command line arguments 
	{
		//we didn't get command line args, so prompt for them
		cout << "Command line arguments error/not provided." << endl;
		cout << "Enter initial altitude/elevation: " << endl;
        cin >> initAlt;
        while (cin.fail())
        {
            cin.clear();
            cin.ignore(); //skip bad input
            cout << "Error - Invalid input." << endl;
            cout << "Enter initial altitude/elevation: " << endl;
            cin >> initAlt;
        }

        cout << "Enter firing angle in degrees (0-90): " << endl;
		cin >> theta;
        while (cin.fail() || theta < 0.0 || theta > 90.0)
        {
            cin.clear();
            cin.ignore(); //skip bad input
            cout << "Error - Invalid input." << endl;
            cout << "Enter firing angle in degrees (0-90): " << endl;
            cin >> theta;
        }

		cout << "Enter initial velocity (m/s): " << endl;
		cin >> initVel;
        while (cin.fail())
        {
            cin.clear();
            cin.ignore(); //skip bad input
            cout << "Error - Invalid input." << endl;
            cout << "Enter initial velocity (m/s): " << endl;
            cin >> initVel;
        }

		cout << "Enter the time step (s) per integration: " << endl;
		cin >> timeStep;
        while (cin.fail() || timeStep < 0.)
        {
            cin.clear();
            cin.ignore(); //skip bad input
            cout << "Error - Invalid input." << endl;
            cout << "Enter the time step (s) per integration: " << endl;
            cin >> timeStep;
        }

		cout << "Enter final time (s): " << endl;
		cin >> duration;
        while (cin.fail() || duration < 0.)
        {
            cin.clear();
            cin.ignore(); //skip bad input
            cout << "Error - Invalid input." << endl;
            cout << "Enter final time (s): " << endl;
            cin >> duration;
        }
	}
	else
	{
	    initAlt = strtod(argv[1], nullptr);
		initVel = strtod(argv[2], nullptr);
		theta = strtod(argv[3], nullptr);
		timeStep = strtod(argv[4], nullptr);
		duration = strtod(argv[5], nullptr);
	}

	cout << "Running Simulation..." << endl;

	// Create filename based on input values
	// Note: stringsteam is a stream for operating on strings. In order
	// to concatinate strings or produce statements like those in java, 
	// we must use this stream. We CANNOT simply concatinate strings 
	// and variables
	stringstream ss;
	ss << "v0_" << initVel << "_ang_" << theta << ".dat";
	string fileName = ss.str();

	ofstream myFile;
	myFile.open(fileName.c_str()); //Must be c style string

	theta *= deg2rad; //firing algle converted to radians

	State projectile;
	projectile.yPos = initAlt;
	projectile.yVel = initVel * sin(theta);
	projectile.yAcc = gravity(projectile.yPos); // g at yPos = 0
	projectile.xPos = 0.0;
	projectile.xVel = initVel * cos(theta);
	projectile.xAcc = a;// a = 0

	if (myFile.fail())
	{
		cout << "Failed To Open File For Writing" << endl;
	}
	else
	{
	    myFile << "InitAlt = " << initAlt << "\tInitVel = " << initVel 
			<< "\ttheta = " << theta << endl;

		double currentTime = 0.0; //current time

		// Do required integration using 4th order Runge Kutta
		while (currentTime <= duration && projectile.yPos >= 0.0)
		{
            currentTime += timeStep;
			rk4(projectile, timeStep);

			if (currentTime <= duration && projectile.yPos >= 0.0)
            {
                // set num of digits after the decimal
                myFile.setf(ios::fixed | ios::showpoint);
                myFile << setprecision(6);
                myFile << currentTime << endl;
                myFile << setprecision(9);
                myFile << "\ty=" << projectile.yPos << " y'=" << 
					projectile.yVel << " y''=" << projectile.yAcc << 
					"\n\tx=" << projectile.xPos << " x'=" << 
					projectile.xVel << " x''=" << projectile.xAcc << 
					endl;
            }
		}
	}

	myFile.close();

	cout << "Simulation complete. Output written to " << fileName << "." << endl;
	return 0;
}

double gravity(const double h)
{
	//return gravity yAcc at current altitude = g * (re/(re+h))^2
	return g * ((re / (re + h)) * (re / (re + h)));
}

void rk4(State & object, const double dt)
{
	//Runge Kutta Method
	//Given an IVP of the form:
	//  y' = f(x, y)
	//  y(xsub[0]) = ysub[0]
	//
	//  xsub[n+1] = xsub[n] + h
	//
	//  ysub[n+1] = ysub[n]+(1/6)(ksub[1]+2(ksub[2]+ksub[3])+ksub[4])
	//  where
	//  ksub[1] = h * f(xsub[n], ysub[n])
	//  ksub[2] = h * f(xsub[n] + h/2, ysub[n] + ksub[1]/2)
	//  ksub[3] = h * f(xsub[n] + h/2, ysub[n] + ksub[2]/2)
	//  ksub[4] = h * f(xsub[n] + h, ysub[n] + ksub[3])

	State k1, k2, k3, k4;

	k1.xPos = object.xPos;
	k1.xVel = object.xVel;
	k1.xAcc = object.xAcc;
	k1.yPos = object.yPos;
	k1.yVel = object.yVel;
	k1.yAcc = object.yAcc;

	k2.xPos = object.xPos + k1.xVel * dt / 2.;
	k2.xVel = object.xVel + k1.xAcc * dt / 2.;
	k2.xAcc = object.xAcc;
	k2.yPos = object.yPos + k1.yVel * dt / 2.;
	k2.yVel = object.yVel + k1.yAcc * dt / 2.;
	k2.yAcc = gravity(k2.yPos);//g recomputed based on yPos

	k3.xPos = object.xPos + k2.xVel * dt / 2.;
	k3.xVel = object.xVel + k2.xAcc * dt / 2.;
	k3.xAcc = object.xAcc;
	k3.yPos = object.yPos + k2.yVel * dt / 2.;
	k3.yVel = object.yVel + k2.yAcc * dt / 2.;
	k3.yAcc = gravity(k3.yPos);//g recomputed based on yPos

	k4.xPos = object.xPos + k3.xVel * dt;
	k4.xVel = object.xVel + k3.xAcc * dt;
	k4.xAcc = object.xAcc;
	k4.yPos = object.yPos + k3.yVel * dt;
	k4.yVel = object.yVel + k3.yAcc * dt;
	k4.yAcc = gravity(k4.yPos);//g recomputed for yPos

	object.xPos += dt*(k1.xVel + 2.*(k2.xVel + k3.xVel) + k4.xVel)/6.;
	object.xVel += dt*(k1.xAcc + 2.*(k2.xAcc + k3.xAcc) + k4.xAcc)/6.;
	object.xAcc += 0;//object.xAcc; //xAcc doesn't change for our sim
	object.yPos += dt*(k1.yVel + 2.*(k2.yVel + k3.yVel) + k4.yVel)/6.;
	object.yVel += dt*(k1.yAcc + 2.*(k2.yAcc + k3.yAcc) + k4.yAcc)/6.;
	object.yAcc = gravity(object.yPos);
}
