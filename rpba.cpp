/*

RPBA - Robust Parallel Bundle Adjustment

File rpba.cpp

Description: Input controlled by file "parameters", call of parallel bundle adjustment and output



Copyright 2019 Helmut Mayer, Bundeswehr University Munich, Germany, Helmut.Mayer@unibw.de

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#include <cstdio>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "rpba.h"
#include "system.h"
#include <iomanip>      // std::setprecision
#include <fstream>      // std::ofstream



using std::string;
using std::ifstream;
using std::ofstream;
using namespace Eigen;



void ComputeRotationMatrixfromRodriguesAngle(const double a, const double b, const double c,
					     PMat &Rabc);
void ComputeRodriguesAnglefromRotationMatrix(double &a, double &b, double &c,
		const PMat &Rabc);


int main(int argc, char** argv) {
	// load parameters from parameter file
	if( argc > 1 ) {
		std::cout << "Parameter file not used\n";
		return 1;
	}


	boost::property_tree::ptree pt;
	boost::property_tree::ini_parser::read_ini("parameters", pt);

	std::string ID = pt.get<std::string>("InputDir");
	std::string File = pt.get<std::string>("File");
	std::string rf = pt.get<std::string>("Robust");
	std::string ee = pt.get<std::string>("Exteval");
	std::string k2i = pt.get<std::string>("k2");
	std::string k4i = pt.get<std::string>("k4");
	std::string fl = pt.get<std::string>("f");
	std::string k24outi = pt.get<std::string>("k24out");
	std::string mti = pt.get<std::string>("Threads");

	int mt = 0;
	bool k2f, k4f, flf, extevalflag, robustflag, k24outflag;

	std::stringstream ees(ee), rfs(rf), k2s(k2i), k4s(k4i), fls(fl), k24outs(k24outi), mts(mti), inputfiles, outputfiles;

	ees >> extevalflag;
	rfs >> robustflag;

	k2s >> k2f;
	k4s >> k4f;
	fls >> flf;

	k24outs >> k24outflag;

	mts >> mt;

	inputfiles << ID << "/" << File;
	std::string inputfile = inputfiles.str();

	std::cout << "Inputfile: " << inputfile << std::endl << std::endl;
	std::cout << "Robustflag: " << robustflag << std::endl << std::endl;

	if (mt < 1)
		mt = System::threadCount();

	std::cout << "Threads: " << mt << std::endl << std::endl;


	std::cout << "K2: " << k2f << std::endl;
	std::cout << "K4: " << k4f << std::endl;
	std::cout << "Focal length: " << flf << std::endl << std::endl;
	std::cout << "K24out: " << k24outflag << std::endl << std::endl;


	std::vector<bool> addflags(7,false);
	addflags[0] = k2f; // k2 - second order coefficient for radial distortion
	addflags[1] = k4f; // k4 - fourth order coefficient for radial distortion
	addflags[2] = flf; // f - focal length or photogrammetric camera constant

	int nraddpar = 0;
	for (int i = 0; i < 3; ++i)
		if (addflags[i])
			++nraddpar;

	MatrixXf XX;
	std::vector< Eigen::MatrixXf > xg;
	std::vector<std::vector<float> > wx, wy, wxy;
	std::vector<std::vector<IMAPNR> > pointXXv;
	std::vector< Eigen::Matrix<float, 3, 4> > PMatrices;
	std::vector<Camera> cameras;
	std::vector<int> iwidth, iheight;

	InputBlock(inputfile, nraddpar, XX, xg,
			wx, wy, wxy,
			pointXXv, PMatrices, cameras,
			iwidth, iheight,
			k24outflag);



	std::ostringstream os;
	const System::Time start = System::getTickCount();

	CalcPLSA(XX, xg,
			wx, wy, wxy,
			pointXXv, PMatrices, cameras,
			iwidth, iheight,
			robustflag,
			addflags,
			os, mt, extevalflag);

	if ((int) PMatrices.size() < 50 || mt == 1)
		std::cout << os.str();

	std::cout << "\n\nBundle adjustment runtime: " << (float) (System::getTickCount()-start) / 1000.f << "s\n\n";



/*	if (!extevalflag) {
		outputfiles << ID << "/out-" << File;
		std::string outputfile = outputfiles.str();

		OutputBlock(outputfile, XX, xg, pointXXv, PMatrices,
				cameras, k24outflag);
	} // if extevalflag */

	
	
/*
	// Output of text file with 3D coordinates
	std::stringstream cloudfiles;
	cloudfiles << ID << "/" << "cloud.txt";
	std::string cloudfile = cloudfiles.str();
	std::cout << "Outputfile: " << cloudfile << std::endl;
	ofstream foutc(cloudfile);
	foutc << std::setprecision(8);
	int counter = 0;
	for (int i = 0; i < XX.cols(); ++i) {
		if (pointXXv[i].size() > 0)
			foutc << XX(0,i) << ' ' << XX(1,i) << ' ' << XX(2,i) << '\n';
		else
			++counter;
	}
	foutc.close();

	// Counter gives number of 3D points which have been deleted during (robust) bundle adjustment
	std::cout << "Counter: " << counter << std::endl;
*/

	return 0;
}



void InputBlock(const std::string &inputfile, const int nraddpar,
		Eigen::MatrixXf &XX,
		std::vector<Eigen::MatrixXf> &xg,
		std::vector<std::vector<float> > &wx, std::vector<std::vector<float> > &wy,
		std::vector<std::vector<float> > &wxy,
		std::vector<std::vector<IMAPNR> > &pointXXv,
		std::vector<PMat> &PMatrices,
		std::vector<Camera> &cameras,
		std::vector<int> &iwidth, std::vector<int> &iheight,
		const bool k24outflag)
{
	ifstream fin(inputfile);

	int nrcameras,nrpoints,nrobs;

	fin >> nrcameras;
	fin >> nrpoints;
	fin >> nrobs;

	std::cout << "Cameras: " << nrcameras << "  Points: " << nrpoints << " Observations " << nrobs << std::endl;

	int nrobs2 = nrobs * 2;
	int nrparimages = (nrcameras - 2) * 6 + 5 + nrcameras * nraddpar;
	int nrunknowns = 3 * nrpoints + nrparimages;
	int redundancy = nrobs2 - nrunknowns;
	std::cout << "Unknowns " << nrunknowns << "  Redundancy " << redundancy << std::endl;

	std::vector<std::vector<float> > xin(nrcameras),yin(nrcameras);

	wx.resize(nrcameras);
	wy.resize(nrcameras);
	wxy.resize(nrcameras);

	pointXXv.resize(nrpoints);

	std::vector<float> minx(nrcameras,1.e20f),miny(nrcameras,1.e20f),maxx(nrcameras,-1.e20f),maxy(nrcameras,-1.e20f);
	int nrcam,nrpoint;
	float x,y;
	IMAPNR point;
	for (int i = 0; i < nrobs; ++i) {
		fin >> nrcam >> nrpoint >> x >> y;
		if (x < minx[nrcam])
			minx[nrcam] = x;
		if (y < miny[nrcam])
			miny[nrcam] = y;
		if (x > maxx[nrcam])
			maxx[nrcam] = x;
		if (y > maxy[nrcam])
			maxy[nrcam] = y;
		xin[nrcam].push_back(x);
		yin[nrcam].push_back(y);
		// set weight to unit matrix
		wx[nrcam].push_back(1.f);
		wy[nrcam].push_back(1.f);
		wxy[nrcam].push_back(0.f);
		point.image = nrcam;
		point.pointnr = (int) xin[nrcam].size() - 1;
		pointXXv[nrpoint].push_back(point);
	}

	PMatrices.resize(nrcameras);
	cameras.resize(nrcameras);
	xg.resize(nrcameras);
	iwidth.resize(nrcameras);
	iheight.resize(nrcameras);

	double a, b, c;
	float f, k2, k4;
	for (int i = 0; i < nrcameras; ++i) {
		cameras[i] = Camera();
		cameras[i].K.setZero();
		fin >> a >> b >> c >> PMatrices[i](0,3) >> PMatrices[i](1,3) >> PMatrices[i](2,3) >> f >> k2 >> k4;
		ComputeRotationMatrixfromRodriguesAngle(a,b,c,PMatrices[i]);
		cameras[i].K(1,1) = cameras[i].K(0,0) = f;
		cameras[i].K(2,2) = 1.f;
		if (k24outflag) {
			cameras[i].distpara[0] = k2;
			cameras[i].distpara[1] = k4;
		} // if
		else {
			cameras[i].distpara[0] = k2 * f * f;
			cameras[i].distpara[1] = k4 * f * f * f * f;
		} // else

		iwidth[i] = 1;
		iheight[i] = 1;

		xg[i].resize(3,(int) xin[i].size());
		for (int j = 0; j < (int) xin[i].size(); ++j) {
			xg[i](0,j) = -xin[i][j];
			xg[i](1,j) = -yin[i][j];
			xg[i](2,j) = 1.f;
		} // for j
	} // for i

	XX.resize(4,nrpoints);
	for (int i = 0; i < nrpoints; ++i) {
		fin >> XX(0,i) >> XX(1,i) >> XX(2,i);
		XX(3,i) = 1.f;
	} // for i

	fin.close();

	return;
}



void OutputBlock(const std::string &outputfile,
		Eigen::MatrixXf &XX,
		std::vector<Eigen::MatrixXf> &xg,
		std::vector<std::vector<IMAPNR> > &pointXXv,
		std::vector<PMat> &PMatrices,
		std::vector<Camera> &cameras,
		const bool k24outflag)
{
	ofstream fout(outputfile);

	int nrobs = 0;

	for (int i = 0; i < (int) pointXXv.size(); ++i)
		nrobs += (int) pointXXv[i].size();

	fout << (int) PMatrices.size() << ' ';
	fout << (int) XX.cols() << ' ';
	fout << nrobs << '\n';

	SBI sortbimage;

	fout << std::setprecision(8);
	for (int i = 0; i < (int) pointXXv.size(); ++i) {
		sort(pointXXv[i].begin(), pointXXv[i].end(), sortbimage);
		for (int j = 0; j < (int) pointXXv[i].size(); ++j) {
			IMAPNR &point = pointXXv[i][j];
			fout << point.image << ' ' << i << ' ' << -xg[point.image](0,point.pointnr) << ' ' << -xg[point.image](1,point.pointnr) << '\n';
		} // for j
	} // for i

	double a, b, c;
	for (int i = 0; i < (int) PMatrices.size(); ++i) {
		ComputeRodriguesAnglefromRotationMatrix(a, b, c, PMatrices[i]);
		fout << a << ' ' << b << ' ' << c << ' ' << PMatrices[i](0,3) << ' ' << PMatrices[i](1,3) << ' ' << PMatrices[i](2,3) << ' ';
		fout << cameras[i].K(0,0) << ' ';
		if (k24outflag)
			fout << cameras[i].distpara[0] << ' ' << cameras[i].distpara[1] << '\n';
		else {
			const double f = cameras[i].K(0,0);
			const double f2 = 1./ (f * f);
			fout << cameras[i].distpara[0] * f2 << ' ' << cameras[i].distpara[1] * f2 * f2 << '\n';
		}
	} // for i

	for (int i = 0; i < XX.cols(); ++i)
		fout << XX(0,i) << ' ' << XX(1,i) << ' ' << XX(2,i) << '\n';

	fout.close();

	return;
} // OutputBlock



// inspired by PBA
void ComputeRotationMatrixfromRodriguesAngle(const double a, const double b, const double c,
		PMat &Rabc){

	if (boost::math::isnan(a) || boost::math::isnan(b) || boost::math::isnan(c))
		std::cout << "NANRabc: " << a << ' ' << b << ' ' << c << '\n';

	const double a2 = a * a, b2 = b * b, c2 = c * c;
	const double ab = a * b, ac = a * c, bc = b * c;

	const double aa = sqrt(a2 + b2 + c2);
	const double ct = aa == 0.?0.5:(1. - cos(aa)) / (aa * aa);
	const double st = aa == 0.?1:sin(aa) / aa;
	Rabc(0,0) = (float) (1. - (b2 + c2) * ct);
	Rabc(0,1) = (float) (ab * ct - c * st);
	Rabc(0,2) = (float) (ac * ct + b * st);
	Rabc(1,0) = (float) (ab * ct + c * st);
	Rabc(1,1) = (float) (1. - (c2 + a2) * ct);
	Rabc(1,2) = (float) (bc * ct - a * st);
	Rabc(2,0) = (float) (ac * ct - b * st);
	Rabc(2,1) = (float) (bc * ct + a * st);
	Rabc(2,2) = (float) (1. - (a2 + b2) * ct);

	return;
} // ComputeRotationMatrixfromRodriguesAngle



// inspired by PBA
void ComputeRodriguesAnglefromRotationMatrix(double &a, double &b, double &c,
		const PMat &Rabc)
{
    const double epsilon0 = 0.01;
    const double epsilon1 = 0.1;
    const double pi = 3.14159265358979323846;

    double trace1 = (Rabc(0,0) + Rabc(1,1) + Rabc(2,2) - 1.) * 0.5;
    if (fabs(Rabc(0,1) - Rabc(1,0)) < epsilon0 &&
    		fabs(Rabc(1,2) - Rabc(2,1)) < epsilon0 &&
		fabs(Rabc(0,2) - Rabc(2,0)) < epsilon0 ) {
    	 	 if (fabs(Rabc(0,1) + Rabc(1,0)) < epsilon1 &&
    	 		fabs(Rabc(1,2) + Rabc(2,1)) < epsilon1 &&
			fabs(Rabc(0,2) + Rabc(2,0)) < epsilon1 && trace1 > 0.9) {
    	 		 a = b = c = 0.;
    	 		 return;
    	 	 }

    	 	 const double xx12 = (Rabc(0,0) + 1.) * 0.5;
    	 	 const double yy12 = (Rabc(1,1) + 1.) * 0.5;
    	 	 const double zz12 = (Rabc(2,2) + 1.) * 0.5;
    	 	 const double xy4 = (Rabc(0,1) + Rabc(1,0)) * 0.25;
    	 	 const double xz4 = (Rabc(0,2) + Rabc(2,0)) * 0.25;
    	 	 const double yz4 = (Rabc(1,2) + Rabc(2,1)) * 0.25;

    	 	 if ((xx12 > yy12) && (xx12 > zz12)) {
    	 		 if (xx12 < epsilon0) {
    	 			 a = 0.;
    	 			 b = c = sqrt(0.5) * pi;
    	 			 return;
    	 		 } // if xx
    	 		 const double t = sqrt(xx12);
    	 		 a = t * pi;
    	 		 b = xy4 / t * pi;
    	 		 c = xz4 / t * pi;
    	 		 return;
    	 	 } // if xx
    	 	 if (yy12 > zz12) {
    	 		 if (yy12 < epsilon0) {
    	 			 a = c = sqrt(0.5) * pi;
    	 			 b = 0.;
    	 			 return;
    	 		 } // if yy
    	 		 const double t = sqrt(yy12);
    	 		 a = xy4 / t * pi;
    	 		 b = t * pi;
    	 		 c = yz4 / t * pi;
    	 		 return;
    	 	 } // if
    	 	 if (zz12 < epsilon0) {
    	 		 a = b = sqrt(0.5) * pi;
    	 		 c = 0.;
    	 		 return;
    	 	 } // if zz
    	 	 const double t  = sqrt(zz12);
    	 	 a = xz4 / t * pi;
    	 	 b = yz4 / t * pi;
    	 	 c = t * pi;
    	 	 return;
    } // if fabs

    const double aa = acos(trace1);
    const double bb = 0.5 * aa / sin(aa);
    a = bb * (Rabc(2,1) - Rabc(1,2));
    b = bb * (Rabc(0,2) - Rabc(2,0));
    c = bb * (Rabc(1,0) - Rabc(0,1));

    return;
} // ComputeRodriguesAnglefromRotationMatrix
