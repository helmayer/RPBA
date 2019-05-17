/*

RPBA - Robust Parallel Bundle Adjustment

File base.h



Copyright 2019 Helmut Mayer, Bundeswehr University Munich, Germany, Helmut.Mayer@unibw.de

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#ifndef _RBASE_H_
#define _RBASE_H_

#include <stdio.h>
#include <vector>
#include <fstream>
#include <ctime>
#include <sys/timeb.h>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <boost/math/special_functions/fpclassify.hpp>

#pragma warning(push)
#pragma warning(disable : 4005)
#include <metis.h>
#pragma warning(pop)


typedef Eigen::Matrix<float,3,4> PMat;		/**< Projection matrix. */
typedef Eigen::Matrix<double,3,4> PMatd;	/**< Projection matrix with double precision. */


#include<Eigen/StdVector>
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(PMat)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(PMatd)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix3d)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix3f)

using std::string;


/**
 * @mainpage RPBA - Robust Parallel Bundle Adjustment
 *
 * @author
 *   Helmut Mayer with contributions by Mario Michelini\n
 *   Visual Computing, Bundeswehr University Munich, Germany
 *
 * @date April 2019
 *
 * @note
 * Please cite
 * \verbatim
@InProceedings{mayer:19,
  author    = "Mayer, H.",
  title     = "{RPBA -- Robust Parallel Bundle Adjustment Based on Covariance Information}",
  pages =        "submitted",
  journal =      "ISPRS Journal of Photogrammetry and Remote Sensing",
  year =         "2019",
} \endverbatim
 */




/**
 * @brief Main function
 *
 * Loads parameter file
 * @param[in] argc     Number arguments
 * @param[in] argv     Argument vector - first element is parameter file
 */
int main(int argc, char** argv);


/** @file */

/** @brief Structure used for pointXXv */
struct IMAPNR{
	/** Image number */
	int image;
	/** Point number */
	int pointnr;
};

/**
 * Camera
 */
class Camera
{
public:
	Eigen::Matrix3f K;				/**< Inner calibration matrix */
	std::vector<float> distpara;	/**< Distortion parameters */

public:
	/**
	 * Creates a new camera
	 */
	Camera() {
		distpara.resize(2, 0);
	}
};

/** @brief Structure for sorting using a custom function object */
struct SBI {
	/**
	 * @brief Function for comparison
	 *
	 * @param[in] i         First element to compare
	 * @param[in] j         Second element to compare
	 *
	 * @result True if preposition holds
	 */
	inline bool operator()(const IMAPNR& i, const IMAPNR& j) const {
		return i.image < j.image;
	}
};



// rpba.cpp
/**
 * @brief Utility function for input of the block
 *
 * @param[in] inputfile            Name of input file
 * @param[in] nraddpar             Number of additional parameters (interior parameters)
 * @param[out] XX                  Homogeneous coordinates of 3D points
 * @param[out] xg                  Homogeneous image coordinates
 * @param[out] wx                  Variance of x-coordinates
 * @param[out] wy                  Variance of y-coordinates
 * @param[out] wxy                 Covariance
 * @param[out] pointXXvio          Links 3D points (row) to pairs of image number and point number (column)
 * @param[out] PMatrices           Projection matrices
 * @param[out] cameras             Cameras (interior parameters)
 * @param[out] iwidth              Vector with image widths
 * @param[out] iheight             Vector with image heights
 * @param[in] k24outflag           Controls if input should be prepared according to the output of PBA
 */
void InputBlock(const std::string &inputfile, const int nraddpar,
		Eigen::MatrixXf &XX,
		std::vector<Eigen::MatrixXf> &xg,
		std::vector<std::vector<float> > &wx, std::vector<std::vector<float> > &wy,
		std::vector<std::vector<float> > &wxy,
		std::vector<std::vector<IMAPNR> > &pointXXv,
		std::vector<PMat> &PMatrices,
		std::vector<Camera> &cameras,
		std::vector<int> &iwidth, std::vector<int> &iheight,
		const bool k24outflag);

/**
 * @brief Utility function for output of the block
 *
 * @param[in] outputfile           Name of output file
 * @param[in] nraddpar             Number of additional parameters (interior parameters)
 * @param[in] XX                   Homogeneous coordinates of 3D points
 * @param[in] xg                   Homogeneous image coordinates
 * @param[in] pointXXvio           Links 3D points (row) to pairs of image number and point number (column)
 * @param[in] PMatrices            Projection matrices
 * @param[in] cameras              Cameras (interior parameters)
 * @param[in] k24outflag           Controls if input should be prepared according to the output of PBA
 */
void OutputBlock(const std::string &outputfile,
		Eigen::MatrixXf &XX,
		std::vector<Eigen::MatrixXf> &xg,
		std::vector<std::vector<IMAPNR> > &pointXXv,
		std::vector<PMat> &PMatrices,
		std::vector<Camera> &cameras,
		const bool k24outflag);





// rbacore.cpp
/** @brief Internal structure for cameras for LSA */
struct LSACamera{
	/** Calibration matrix */
	Eigen::Matrix3d K;
	/** Distortion parameters */
	std::vector<double> distpara;
};



/** @brief Robust bundle adjustment
 * @author Helmut Mayer
 */
class LSA{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	/**
	 * @brief Default constructor
 	 */
	/**
	 * @brief Constructor
	 *
	 * @param[in] cameras			   Camera set
	 * @param[in] addflags             Controls which additional parameters are to be used
	 * @param[in] imgindflagin         Controls if a selected number of images as given in imgindin is used
	 * @param[out] imgindin            Selected number of images for bundle adjustment
	 * @param[in,out] PMatricesi       Projection matrices
	 * @param[in,out] pointXXvio       Links 3D points (row) to pairs of image number and point number (column)
	 * @param[in,out] XXi              Homogeneous coordinates of 3D points
	 * @param[in] XXnrin               Number of 3D points - needed as in the two and three image case XX is of constant length for reasons of efficiency
	 * @param[in] xgin                 Homogeneous image coordinates
	 * @param[in] wx                   Variance of x-coordinates
	 * @param[in] wy                   Variance of y-coordinates
	 * @param[in] wxy                  Covariance
	 * @param[in] W3iji                3D covariance matrix for GCPs
	 * @param[in] gcpflagin            Controls if ground control point (GCP) information is to be used
	 * @param[in] gcpinfo            Vector controlling for which points GCP information is to be used
	 * @param[in] gcptypein            GCP type - 0: xyz, 1: xy, 2: z
	 * @param[in] gcpdatain            GCP data (xyz)
	 * @param[in] ihw2i                Half average image width and height
	 * @param[in] robustflagi          Controls if robust re-weighting is to be used
 	 */
	LSA(const std::vector<Camera> &cameras,
			const std::vector<bool> &addflags,
			const bool imgindflagin, const std::vector<int> &imgindin,
			std::vector<PMatd> &PMatricesi,
			std::vector<std::vector<IMAPNR> > &pointXXvio,
			Eigen::MatrixXd &XXi, const int XXnrin,
			const std::vector<Eigen::MatrixXf> &xgin,
			const std::vector<std::vector<float> > &wx, const std::vector<std::vector<float> > &wy, const std::vector<std::vector<float> > &wxy,
			const std::vector<std::vector<float> > &W3iji,
			const bool gcpflagin, const std::vector<int> &gcpinfo, const std::vector<int> &gcptypein, const std::vector<std::vector<double> > &gcpdatain,
			const float ihw2i,
			const bool robustflagi);


	/** @brief Internal structure in LSA for pointXXv
	 *
	 * Accounts for that indices have to be local */
	struct IMPNR{
		/** Image number (local index) */
		int iimage;
		/** Point number */
		int pointnr;
	};


	/**
	 * @brief Input function for parallel bundle adjustment
	 *
	 * @param[in] addflagsin           Controls which additional parameters are to be used
	 * @param[in] cameras			   Camera set
	 * @param[out] imgindin            Selected number of images for bundle adjustment
	 * @param[in,out] pointXXvio       Links 3D points (row) to pairs of image number and point number (column)
	 * @param[in] PMatricesi           Projection matrices
	 * @param[in] XXi                  Homogeneous coordinates of 3D points
	 * @param[in] XXnri                Number 3D points
	 * @param[in] gcpinfo            Vector controlling for which points GCP information is to be used
 	 */
	void InputLSAP(const std::vector<bool> &addflagsin,
			const std::vector<Camera> &cameras,
			const std::vector<int> &imgindin, std::vector<std::vector<IMAPNR> > &pointXXvio,
			const std::vector<PMatd> &PMatricesi,
			const Eigen::MatrixXd &XXi, const int XXnri,
			const std::vector<int> &gcpinfo);

	/**
	 * @brief Adjustment (where the work is done)
	 *
	 * @param[in] robustflagi          Controls if robust re-weighting is to be used
	 * @param[out] os                  Output string containing logging information
	 * @param[in] itbreakflagi         Controls if compuation is stoped after one iteration of outer loop
	 * @param[in] itmaxnumin           Maximum number of iterations
	 * @param[in] W3iji                3D covariance matrices for GCPs
	 * @param[in] gcpflagin            Controls if ground control point (GCP) information is to be used
	 * @param[in] gcptypeit            GCP type - 0: xyz, 1: xy, 2: z
	 * @param[in] gcpdatait            GCP data (xyz)
	 */
	void Adjust(const bool robustflagi,
			std::ostream& os, const bool itbreakflagi, const int itmaxnumin, const std::vector<std::vector<float> > &W3iji,
			const bool gcpflagin, const std::vector<int> &gcptypeit, const std::vector<std::vector<double> > &gcpdatait);
	/**
	 * @brief Global output function for parallel processing
	 *
	 * @param[out] PMatricesi        Projection matrices
	 * @param[out] pointXXvio       Links 3D points (row) to pairs of image number and point number (column)
	 * @param[out] XXio               Homogeneous coordinates of 3D points
	 * @param[out] gcpinfoo            Vector controlling for which points GCP information is to be used
	 * @param[out] cameraso           Calibration information per camera (not image)
 	 */
	void OutputLSAP(std::vector<PMatd> &PMatricesi,
			std::vector<std::vector<IMAPNR> > &pointXXvio, Eigen::MatrixXd &XXio,
			std::vector<int> &gcpinfoo, std::vector<Camera> &cameraso);

	/**
	 * @brief Global output function
	 *
	 * @param[out] PMatricesi        Projection matrices
	 * @param[out] pointXXvio       Links 3D points (row) to pairs of image number and point number (column)
	 * @param[out] XXio               Homogeneous coordinates of 3D points
 	 * @param[out] cameraso            Calibration information per camera (not image)
	 * @return						  Success
 	 */
	/**  @{ */
	bool OutputLSAc(std::vector<PMat> &PMatricesi,
			std::vector<std::vector<IMAPNR> > &pointXXvio, Eigen::MatrixXf &XXio,
			std::vector<Camera> &cameraso);


	/**
	 * @brief Preparation of cameras
 	 */
	void PrepareCameras();
	/**
	 * @brief Preparation for iteration
	 *
	 * @param[in] cset				  Camera set
	 * @param[in] addflagsin             Controls which additional parameters are to be used
	 * @param[in] gcptypein           GCP type - 0: xyz, 1: xy, 2: z
	 * @param[in] gcpdatain           GCP data (xyz)
 	 */
	void PrepareIteration(const std::vector<Camera> &cset,
			const std::vector<bool> &addflagsin,
			const std::vector<int> &gcptypein, const std::vector<std::vector<double> > &gcpdatain);
	/**
	 * @brief Prepares weights for 2D observations
	 *
	 * @param[in] it            Number of iteration
	 * @param[in] robustflagi    Controls if robust re-weighting is to be used
 	 */
	void PrepareWeights2D(const int it, const bool robustflagi);
	/**
	 * @brief Prepares weights for 3D observations
	 *
	 * @param[in] robustflagi    Controls if robust re-weighting is to be used
	 * @param[in] itbreakflagi         Controls if compuation is stoped after one iteration of outer loop
	 * @param[in] W3ijit               3D covariance matrices for GCPs
 	 */
	void PrepareWeights3D(const bool robustflagi, const bool itbreakflagi, const std::vector<std::vector<float> > &W3ijit);
	/**
	 * @brief Reweights the 2D observations
	 *
	 * 	 * @param[in] iflag		Flag for convergence
 	 */
	void Reweight2D(const int iflag);
	/**
	 * @brief Encodes camera information as unknowns (e.g., transformation or rotation matrix to quaternions)
 	 */
	void EncodeCameras();
	/**
	 * @brief Store results for improved average standard deviation
 	 */
	void StoreResultsIteration();
	/**
	 * @brief Restore results after break
	 *
 	 */
	void RobustDelete();
	void CalcCorrections(const bool solveflag, const double multfact1, const std::vector<double> &XM, const bool srmvflag,
			const bool derivflag, const bool robustflagi,
			const bool itbreakflagi, const std::vector<std::vector<float> > &W3ijit,
			const std::vector<int> &gcptypeit, const std::vector<std::vector<double> > &gcpdatait);
	/** @brief Computation of derivatives of the observations
	 *
	 * @param[in] solveflag             If true, solution vector is computed
	 * @param[in] camerak               Index of camera
	 * @param[in] k                     Index
	 * @param[in] ks                    Local index for images
	 * @param[in] nro                   Number observation (global)
	 * @param[in] xs3                   Inverted third component of a point
	 * @param[in] diff1                 First correction
	 * @param[in] diff2                 Second correction
	 */
	void CalcDerivatives(const bool solveflag, const int camerak,
			const int k, const int ks, const int nro, const double xs3,
			const double diff1, const double diff2);

	/** @brief Computation of calibrated derivatives of the observations
	 *
	 * @param[in,out] A                matrix with derivatives
	 * @param[in] camerak              Index of camera
	 */
	void eCalcCalibratedDerivativesX(Eigen::MatrixXd &A, const int camerak);

	/** @brief Computation of calibrated derivatives of the observations
	 *
	 * @param[in,out] A                matrix with derivatives
	 * @param[in] camerak              Index of camera
	 */
	void eCalcCalibratedDerivatives3(Eigen::Matrix<double, 2, 3> &A, const int camerak);


	/** @brief Computation of 3D derivatives for GCPs
	 *
	 * @param[in] i                     Index of point i
	 * @param[in] gcpnum                Index of GCP
	 * @param[in] robustflagi           Controls if robust re-weighting is to be used
	 * @param[in] itbreakflagi          Controls if compuation is stoped after one iteration of outer loop
	 * @param[in] W3ijit                3D covariance matrices for GCPs
	 * @param[in] gcptypeit             GCP type - 0: xyz, 1: xy, 2: z
	 * @param[in] gcpdatait             GCP data (xyz)
	 * @param[in] derivflag             If true, derivatives are computed
	 */
	void CalcDerivatives3DGCP(const int i, const int gcpnum,
			const bool robustflagi, const bool itbreakflagi, const std::vector<std::vector<float> > &W3ijit,
			const std::vector<int> &gcptypeit, const std::vector<std::vector<double> > &gcpdatait, const bool derivflag);


	/** @brief Decoding of cameras given the unknowns
	 *
	 * @param[in] XM                    Vector with unknowns of bundle adjustment
	 */
	void DecodeCameras(const std::vector<double> &XM);

	/** @brief Sparse Conjugate Gradients
	 *
	 * @param[in] begini                Controls with which image to start - for absolute orientation solution this is zero, for image only solution it is one
	 * @return			   		True if error
	 */
	bool SparsePCG(const int begini);
	/** @brief Core of sparse bundle adjustment
	 *
	 * Bundle solution separated for 3D points, projection matrices and additional parameters. The algorithm follows the book
	 * of Mikhail et al. "Introduction to Modern Photogrammetry".
	 *
	 * @param[in] multfact              Multiplication factor for main diagonal of normal equation system according to Levenberg-Marquardt
	 * @param[in] robustflagi           Controls if robust re-weighting is to be used
	 * @param[in] itbreakflagi         Controls if compuation is stoped after one iteration of outer loop
	 * @param[in] W3ijit               3D covariance matrices for GCPs
	 * @param[in] gcptypeit            GCP type - 0: xyz, 1: xy, 2: z
	 * @param[in] gcpdatait            GCP data (xyz)
	 * @return					True if error
	 */
	bool CalcBundleSep(const double multfact, const bool robustflagi,
			const bool itbreakflagi, const std::vector<std::vector<float> > &W3ijit,
			const std::vector<int> &gcptypeit, const std::vector<std::vector<double> > &gcpdatait);
	/** @brief Converts (external) pointXXvio to internal pointXXv
	 *
	 * @param[in] XXi              Homogeneous coordinates of 3D points
	 * @param[in] pointXXvi            Rows consist of 3D points and columns of pairs of point and (global) image numbers
	 * @param[in] gcpinfoin            Vector controlling for which points GCP information is to be used
	 */
	void PointXXvIn(const Eigen::MatrixXd &XXi, std::vector<std::vector<IMAPNR> > &pointXXvi, const std::vector<int> &gcpinfoin);
	/** @brief Converts internal pointXXv after bundle adjustment to pointXXvio
	 *
	 * @param[out] XXio               Homogeneous coordinates of 3D points
	 * @param[out] pointXXvio         Rows consist of 3D points and columns of pairs of point and (global) image numbers - global version
	 * @param[out] gcpinfoo            Vector controlling for which points GCP information is to be used
	 */
	void PointXXvOut(Eigen::MatrixXd &XXio, std::vector<std::vector<IMAPNR> > &pointXXvio, std::vector<int> &gcpinfoo);
	/** @brief Clears memory at the end of LSA */


	protected:
	/** True for initialization */
	bool iniflag;
	/** True if index in imgind is used */
	bool imgindflag;
	/** Flag for longer focal length */
	bool lflflag;


	/** Number of points */
	int imagepointnr;
	/** Number of 3D points */
	int nrpts;
	/** Number of images */
	int nrimages;
	/** Number of cameras */
	int nrcameras;
	/** Number of relevant points in XX */
	int XXnr;

	/** Offset values @{ */
	int nrioff,nrioff3; /**@}*/

	/** Number of parameters per image (reduced) @{ */
	int nrpari,nrparir; /**@}*/

	/** Number additional parameters @{ */
	int nraddpi; /**@}*/
	/** Vector with flags controlling which additional parameters to estimate */
	std::vector<bool> addflags;

	/** Minimum number of images a point needs to be seen in */
	int minimagenr;

	/** Minimum number of images a point needs to be seen in after deletion during robust estimation */
	int minimagenrdel;

	/** Number of observations (2 * imagepointnr) @{ */
	int nrobs, nrobsnew; /**@}*/

	/** Number of unknowns @{ */
	int nrunknowns, nrunknownsnew; /**@}*/

	/** Redundancy @{ */
	int redundancy, redundancynew; /**@}*/

	/** Maximum number of image observations for a 3D point - used to scale the amount of memory needed */
	int maxnrobsc;

	/**         Controls with which image to start - for absolute orientation solution this is zero, for image only solution it is one */
	int begini;

	/** Index of camera for which scale is kept constant - usually, maximum distance of all cameras to first camera is employed */
	int normcamera;
	/** Index to coordinate with maximum value for normcamera */
	int esvmaxm;
	/** Indices of coordinates for normcamera */
	std::vector<int> esvm;
	/** Half average image width and height */
	const float ihw2;
	/** Threshold for wrong projections */
	const double defdiff;


	/** Controls number of standard deviations for robust estimation */
	const float kf = 9.f;

	/** Global sum of weights - used to compute normalized sigma0 values */
	float wsumg;

	/** Sum of weights - used to compute normalized sigma0 values */
	float wsum;

	/** Sum of weights - used to compute normalized sigma0 values for absolute orientation */
	float wsumao;

	/** Average standard deviation for original weights @{ */
	float sig0, sig0old, sig0rw, sig0in, sig0ao, sig0aoi, sig0aoiold; /**@}*/

	/** Vector of unknowns @{ */
	std::vector<double> XPMatrices, XPMatricesnew; /**@}*/

	/** Derived image coordinates (homogeneous) */
	std::vector<double> x;
	/** Homogeneous 3D coordinates for one point */
	Eigen::Vector4d eXX4;
	/** Homogeneous 3D coordinates */
	Eigen::MatrixXd XX;
	/** Homogeneous 3D coordinates */
	Eigen::MatrixXd XXnew;

	/** Matrix for projection of 3D point @{ */
	Eigen::Matrix3d eR;
	Eigen::Vector3d ex3, ex3c, ex3u; /**@}*/


	/** Vector of LSAcamera @{ */
	std::vector<LSACamera> cameras,camerasold,camerasnew,camerasgiven; /**@}*/

	/** Flags controlling if elements of normal equation matrix are used @{ */
	std::vector<std::vector<bool> > Npqqnumv; /**@}*/

	/** True if absolute orientation, i.e., 6 parameters for all cameras, is to be determined */
	bool aoflag;
	/** Number of 3D observations */
	int nrobsao;
	/** Number of 3D points (GCP) */
	int nrao;

	/** True if information for ground control points (GCP) is available */
	bool gcpflag;
	/** Vector controlling for which points GCP information is to be used */
	std::vector<int> gcpinfo;
	/** GCP type - 0: xyz, 1: xy, 2: z */
	std::vector<int> gcptype;
	/** GCP data */
	std::vector<std::vector<double> > gcpdata;
	/** Number of observations for GCPs */
	int nrobsgcp;
	/** Number of GCPs */
	int nrgcp;

	/** Vector with indices to camera parameters */
	std::vector<int> camparv;

	/** Index from valid points to 3D points in XX */
	std::vector<int> nrpv;

	/** Number observations per image summed up (old) @{ */
	std::vector<int> nrobv, nrobvo; /**@}*/

	/** Relates images to LSAcamera */
	std::vector<int> cameraind;

	/** Vector with indices to images, to restrict adjustment to certain images */
	std::vector<int> imgind;

	/** Indices used in CalcBundleSep @{ */
	std::vector<int> nbeginv,nendv,colv; /**@}*/

	/** Weights as inverse variances and covariances @{ */
	const std::vector<std::vector<float> > wxi, wyi, wxyi; /**@}*/

	/** Weight matrix */
	std::vector<std::vector<float > > Wij;

	/** Weight matrix 3D @{ */
	std::vector<std::vector<float> > W3ij,W3ijold; /**@}*/

	/** Weight matrix @{ */
	std::vector<std::vector<std::vector<float> > > Wijoldv; /**@}*/

	/** Flag if observation is to be deleted @{ */
	std::vector<bool> ISflag, ISflagnew; /**@}*/

	/** Vector used to store points to be erased @{ */
	std::vector<int> erasev, erasevpts; /**@}*/

	/** Vector with standard deviation per image @{ */
	std::vector<float> sig0rmv,sig00v;
	std::vector<std::vector<float> > sig0rmvvs; /**@}*/

	/** Vector with standard deviation per 3D points @{ */
	std::vector<float> sig0aorv,sig0aorvs; /**@}*/

	/** Projection matrix */
	std::vector<PMatd> PMatrices;

	/** Rotation matrix */
	std::vector<Eigen::Matrix3d> eRM;
	/** New rotation matrix */
	std::vector<Eigen::Matrix3d> eRMnew;
	/** Projection center */
	Eigen::MatrixXd ePC;
	/** New projection center */
	Eigen::MatrixXd ePCnew;

	/** Measured image coordinates */
	const std::vector<Eigen::MatrixXf> &xg;

	/** Matrices used for the computation of corrections @{ */
	Eigen::Matrix3d eRabc, eRMold0;
	Eigen::Vector3d eesd, eSResd, ePCold0;
	/**@}*/

	/** 3D matrices  @{ */
	Eigen::Matrix3d eNppi, eNppi1;
	/**@}*/
	/** 3D vectors */
	Eigen::Vector3d etppi;

	/** Matrices used in core of sparse bundle adjustment @{ */
	Eigen::MatrixXd eNlij3, eBpi2, eBpi2r, eBpi2t, eBpi2rt, eBpi3, eBpi3r, eBpi3t, eBpi3rt;
	Eigen::MatrixXd eBpi23, eBpit23, efp2;
	Eigen::MatrixXd eBhij, eBhijt, eBhij2, eBhij2t, eNlijN;
	/**@}*/

	/** Vector(s) of Matrix used in core of sparse bundle adjustment @{ */
	std::vector<std::vector<Eigen::MatrixXd> > eNpqqv, eNpqqvni, eNpqqvnew;

	std::vector<Eigen::MatrixXd> eNlij;
	std::vector<Eigen::MatrixXd> exiv, exivold, eCv, erkv, erk1v, edkv, eAdkv, ehkv, ehk1v, ehkdv; /**@}*/

	/** Links 3D points (row) to pairs of image number and point number IMPNR (column) */
	std::vector<std::vector<IMPNR> > pointXXv;

}; //LSA




// rpbacore.cpp
/**
 * @brief Utility function for parallel bundle adjustment
 *
 * @param[in,out] XX               Homogeneous coordinates of 3D points
 * @param[in] xg                   Homogeneous image coordinates
 * @param[in] wx                   Variance of x-coordinates
 * @param[in] wy                   Variance of y-coordinates
 * @param[in] wxy                  Covariance
 * @param[in,out] pointXXvio       Links 3D points (row) to pairs of image number and point number (column)
 * @param[in,out] PMatrices        Projection matrices
 * @param[in,out] cameras          Cameras
 * @param[in] iwidth               Vector with image widths
 * @param[in] iheight              Vector with image heights
 * @param[in] robustflag           Controls if robust re-weighting is to be used
 * @param[in] addflagsin           Controls which additional parameters are to be used
 * @param[out] os                  Output stream
 * @param[in] mt                   Number of threads
 * @param[in] extevalflag	    evaluation of external block
 */
void CalcPLSA(Eigen::MatrixXf &XX,
		std::vector<Eigen::MatrixXf> &xg,
		const std::vector<std::vector<float> > &wx, const std::vector<std::vector<float> > &wy,
		const std::vector<std::vector<float> > &wxy,
		std::vector<std::vector<IMAPNR> > &pointXXvio,
		std::vector<PMat> &PMatrices,
		std::vector<Camera> &cameras,
		const std::vector<int> &iwidth, const std::vector<int> &iheight,
		const bool robustflag,
		const std::vector<bool> &addflagsin,
		std::ostream& os, const int mt, const bool extevalflag);

/**
 * @brief Compute corrections for 3D point for multiple images
 *
 * @param[in] cameraps             Cameras
 * @param[in] X0                   X-coordinate of 3D point
 * @param[in] X1                   Y-coordinate of 3D point
 * @param[in] X2                   Z-coordinate of 3D point
 * @param[in] pxv                  Links 3D points (row) to pairs of image number and point number (column)
 * @param[out] corrv               Vector of corrections used in adjustment
 * @param[in] PM                   Projection matrices
 * @param[in] xg                   Homogeneous image coordinates
 * @param[in] wx                   Variance of x-coordinates
 * @param[in] wy                   Variance of y-coordinates
 * @param[in] wxy                  Covariance
 * @param[in] wxo                  Original variance of x-coordinates
 * @param[in] wyo                  Original variance of y-coordinates
 * @param[in] wxyo                 Original covariance
 * @param[out] sig0                Standard deviation
 * @param[out] sig0v               Vector of squared errors per image
 * @param[out] p1v                 Vector of derivatives of first row of PM
 * @param[out] p2v                 Vector of derivatives of second row of PM
 * @param[out] p3v                 Vector of derivatives of third row of PM
 */
void CalcCorrPointMulti(const std::vector<Camera> &cameraps,
		const double X0, const double X1, const double X2,
		std::vector<IMAPNR> &pxv,
		std::vector<double> &corrv,
		const std::vector<PMatd> &PM,
		const std::vector<Eigen::MatrixXf> &xg,
		const std::vector<float> &wx, const std::vector<float> &wy, const std::vector<float> &wxy,
		const std::vector<float> &wxo, const std::vector<float> &wyo, const std::vector<float> &wxyo,
		float &sig0, std::vector<float> &sig0v,
		std::vector<float> &p1v, std::vector<float> &p2v, std::vector<float> &p3v);

/**
 * @brief Calibrate derivatives
 *
 * @param[in, out] xp0             Bpp2[0][0]
 * @param[in, out] yp0             Bpp2[1][0]
 * @param[in, out] xp1             Bpp2[0][1]
 * @param[in, out] yp1             Bpp2[1][1]
 * @param[in, out] xp2             Bpp2[0][2]
 * @param[in, out] yp2             Bpp2[1][2]
 * @param[in] x                    x-coordinate
 * @param[in] y                    y-coordinate
 * @param[in] k2                   Radial distortion square term
 * @param[in] k4                   Radial distortion quartic term
 * @param[in] fx                   Camera constant in x
 * @param[in] fy                   Camera constant in y
 * @param[in] s                    Skew
 */
void CalcCalibratedDerivatives(double &xp0, double &yp0, double &xp1, double &yp1, double &xp2, double &yp2,
		const double x, const double y, const double k2, const double k4,
		const double fx, const double fy, const double s);

/**
 * @brief Bundle adjustment for individual point - 3D intersection
 *
 * @param[in] cameraps          	  Cameras
 * @param[in] pi	          	  	  Index of point in XX / pointXXv
 * @param[in,out] XX               Homogeneous coordinates of 3D points
 * @param[in] xg                   Homogeneous image coordinates
 * @param[in] wx                   Variance of x-coordinates
 * @param[in] wy                   Variance of y-coordinates
 * @param[in] wxy                  Covariance
 * @param[in,out] pointXXvio       Links 3D points (row) to pairs of image number and point number (column)
 * @param[in,out] PMatrices        Projection matrices
 * @param[in] imgindbv             Inverse index for image number
 * @param[in] nrthreads            Number threads
 * @param[in] pitflags             Flag is true if there is a camera in the corresponding thread
 * @param[out] wXXv                variance in X-direction of 3D point per thread
 * @param[out] wYYv                variance in Y-direction of 3D point per thread
 * @param[out] wZZv                variance in Z-direction of 3D point per thread
 * @param[out] wXYv                covariance in X-Y-direction of 3D point per thread
 * @param[out] wXZv                covariance in X-Z-direction of 3D point per thread
 * @param[out] wYZv                covariance in Y-Z-direction of 3D point per thread
 * @param[out] sigio               Standard deviation
 * @param[out] redundancy          Redundancy
 * @param[out] sig0vv              Vector of standard deviations per image
 * @param[out] sig0wv              Vector of standard deviations per image derived per median and used in robust estimation
 * @param[in] robustmode           Controls different modes of robust estimation
 * @param[in] weightonlyflag       If true, only weights are computed (and no improved XX values)
 * @param[in] wxyflag	          Controls if 3D covariance matrix per thread is computed
 */
bool CalcLSAXMulti(const std::vector<Camera> &cameraps,
		const int pi, Eigen::MatrixXd &XX,
		const std::vector<Eigen::MatrixXf> &xg,
		const std::vector<std::vector<float> > &wx, const std::vector<std::vector<float> > &wy,
		const std::vector<std::vector<float> > &wxy,
		std::vector<std::vector<IMAPNR> > &pointXXvio,
		const std::vector<PMatd> &PMatrices,
		const std::vector<int> &imgindbv, const int nrthreads, const std::vector<bool> &pitflags,
		std::vector<float> &wXXv, std::vector<float> &wYYv, std::vector<float> &wZZv,
		std::vector<float> &wXYv, std::vector<float> &wXZv, std::vector<float> &wYZv,
		double &sig0io, long &redundancy, std::vector<std::vector<float> > &sig0vv,
		const std::vector<float> &sig0wv, const int robustmode, const bool weightonlyflag, const bool wxyflag);

/**
 * @brief Evaluation function computing back projection error for an image block
 *
 * @param[in] nraddp          	  Number of additional parameters used for computing the redundancy
 * @param[in] cameraps          	  Cameras
 * @param[in] XXd                  Homogeneous coordinates of 3D points
 * @param[in] xg                   Homogeneous image coordinates
 * @param[in] wx                   Variance of x-coordinates
 * @param[in] wy                   Variance of y-coordinates
 * @param[in] wxy                  Covariance
 * @param[in] pointXXvio           Links 3D points (row) to pairs of image number and point number (column)
 * @param[in] PMatricesd           Projection matrices
 * @param[in] ihw2                 Scaling factor for standard deviation consisting of the average width and height of an image
 */
void EvalLSA(const int nraddp, const std::vector<Camera> &cameraps,
		const Eigen::MatrixXd &XXd,
		const std::vector<Eigen::MatrixXf> &xg,
		const std::vector<std::vector<float> > &wx, const std::vector<std::vector<float> > &wy,
		const std::vector<std::vector<float> > &wxy,
		std::vector<std::vector<IMAPNR> > &pointXXv,
		const std::vector<PMatd> &PMatricesd,
		const float ihw2);

/**
 * @brief Function to partition an image block into sub-blocks based on METIS
 *
 * @param[in] pointXXv             Links 3D points (row) to pairs of image number and point number (column)
 * @param[in] nrimages             Number of images in blocks
 * @param[in] nrthreads            Number of sub-blocks the block is to be splitted
 * @param[out] part                Index to which sub-block each image belongs to
 * @param[out] imagesperpart       Number of images per sub-block
 */
void PartitionpointXXv(std::vector<std::vector<IMAPNR> > &pointXXv, const int nrimages,
		const int nrthreads, std::vector<idx_t> &part,
		std::vector<std::vector<int> > &imagesperpart);





// Mathematical functions
/**
 * @brief Applies radial distortion to projected image coordinates
 *
 * @param[in,out] x                   x-coordinate of point
 * @param[in,out] y                   y-coordinate of point
 * @param[in] camera                  Pointer to camera
 */
template<class T, class C>
inline void RadialReconstructionImage(T& x, T& y, const C& camera)
{
	const double dist = x * x + y * y;
	const double z = 1. +  (T)(camera.distpara[0] * dist + camera.distpara[1] * dist * dist);
	x *= (T) z;
	y *= (T) z;
}

/**
 * @brief Rotation matrix from Rodrigues parameters
 *
 * @param[in] a                   First component of rotation
 * @param[in] b                   Second component of rotation
 * @param[in] c                   Third component of rotation
 * @param[out] Rabc               Rotation matrix
 */
void ComputeRotationMatrixfromRodrigues(const double a, const double b, const double c,
		Eigen::Matrix3d &Rabc);

/**
 * @brief Rectify rotation matrix based on Quaternions
 *
 * @param[out] R               Rotation matrix
 */
// Rectify rotation matrix based on quaternions
template<class T>
void RectifyRotationMatrix3(Eigen::Matrix<T, 3, 3> &R) {

  T qv[4],a[3];

  // We use quaternions for the representation. This avoids
  // mostly trigonometric functions and is independent of
  // definitions of axes and angles
  // Foerstner paper Ebner Festschrift 1999 p. 86
  a[0] = R(2,1) - R(1,2);
  a[1] = R(0,2) - R(2,0);
  a[2] = R(1,0) - R(0,1);
  if (fabs(a[0]) < 1.e-20)
	  a[0] = 1.e-20;
  if (fabs(a[1]) < 1.e-20)
	  a[1] = 1.e-20;
  if (fabs(a[2]) < 1.e-20)
	  a[2] = 1.e-20;
  const T al = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
  a[0] /= al; a[1] /= al; a[2] /= al;
  const T phi = atan2(al,R(0,0) + R(1,1) + R(2,2) - 1.);
  // vector of quaternions
  qv[0] = cos(phi * 0.5);
  for (int i = 1; i < 4; ++i)
    qv[i] = a[i-1] * sin(phi * 0.5);

  R(0,0) = qv[0] * qv[0] + qv[1] * qv[1] - qv[2] * qv[2] - qv[3] * qv[3];
  R(0,1) = 2 * (qv[1] * qv[2] - qv[0] * qv[3]);
  R(0,2) = 2 * (qv[1] * qv[3] + qv[0] * qv[2]);
  R(1,0) = 2 * (qv[2] * qv[1] + qv[0] * qv[3]);
  R(1,1) = qv[0] * qv[0] - qv[1] * qv[1] + qv[2] * qv[2] - qv[3] * qv[3];
  R(1,2) = 2 * (qv[2] * qv[3] - qv[0] * qv[1]);
  R(2,0) = 2 * (qv[3] * qv[1] - qv[0] * qv[2]);
  R(2,1) = 2 * (qv[3] * qv[2] + qv[0] * qv[1]);
  R(2,2) = qv[0] * qv[0] - qv[1] * qv[1] - qv[2] * qv[2] + qv[3] * qv[3];

  return;
} // RectifyRotationMatrix



/**
 * @brief Rectify rotation matrix based on Quaternions
 *
 * @param[out] R               Rotation matrix
 */
// Rectify rotation matrix based on quaternions
template<class T>
void RectifyRotationMatrix4(Eigen::Matrix<T, 3, 4> &R) {

  T qv[4],a[3];

  // We use quaternions for the representation. This avoids
  // mostly trigonometric functions and is independent of
  // definitions of axes and angles
  // Foerstner paper Ebner Festschrift 1999 p. 86
  a[0] = R(2,1) - R(1,2);
  a[1] = R(0,2) - R(2,0);
  a[2] = R(1,0) - R(0,1);
  if (fabs(a[0]) < 1.e-20)
	  a[0] = 1.e-20;
  if (fabs(a[1]) < 1.e-20)
	  a[1] = 1.e-20;
  if (fabs(a[2]) < 1.e-20)
	  a[2] = 1.e-20;
  const T al = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
  a[0] /= al; a[1] /= al; a[2] /= al;
  const T phi = atan2(al,R(0,0) + R(1,1) + R(2,2) - 1.);
  // vector of quaternions
  qv[0] = cos(phi * 0.5);
  for (int i = 1; i < 4; ++i)
    qv[i] = a[i-1] * sin(phi * 0.5);

  R(0,0) = qv[0] * qv[0] + qv[1] * qv[1] - qv[2] * qv[2] - qv[3] * qv[3];
  R(0,1) = 2 * (qv[1] * qv[2] - qv[0] * qv[3]);
  R(0,2) = 2 * (qv[1] * qv[3] + qv[0] * qv[2]);
  R(1,0) = 2 * (qv[2] * qv[1] + qv[0] * qv[3]);
  R(1,1) = qv[0] * qv[0] - qv[1] * qv[1] + qv[2] * qv[2] - qv[3] * qv[3];
  R(1,2) = 2 * (qv[2] * qv[3] - qv[0] * qv[1]);
  R(2,0) = 2 * (qv[3] * qv[1] - qv[0] * qv[2]);
  R(2,1) = 2 * (qv[3] * qv[2] + qv[0] * qv[1]);
  R(2,2) = qv[0] * qv[0] - qv[1] * qv[1] - qv[2] * qv[2] + qv[3] * qv[3];

  return;
} // eRectifyRotationMatrix






#endif
// _RBASE_H_
