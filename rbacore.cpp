/*

RPBA - Robust Parallel Bundle Adjustment

File rbacore.cpp

Description: Core bundle adjustment functionality



Copyright 2019 Helmut Mayer, Bundeswehr University Munich, Germany, Helmut.Mayer@unibw.de

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/


#include "rpba.h"



inline void eAAdv(const int begini, const int nb,
		const std::vector<std::vector<Eigen::MatrixXd> > &A,
		const std::vector<std::vector<bool> > &Anum,
		const std::vector<Eigen::MatrixXd> &b,
		std::vector<Eigen::MatrixXd> &c);
inline void evvdv(const int begini, const int nb, const std::vector<Eigen::MatrixXd> &a,
		const std::vector<Eigen::MatrixXd> &b,
		double &c);



LSA::LSA(const std::vector<Camera> &cameras,
		const std::vector<bool> &addflagsin,
		const bool imgindflagin, const std::vector<int> &imgindin,
		std::vector<PMatd> &PMatricesi,
		std::vector<std::vector<IMAPNR> > &pointXXvio,
		Eigen::MatrixXd &XXi, const int XXnrin,
		const std::vector<Eigen::MatrixXf> &xgin,
		const std::vector<std::vector<float> > &wx, const std::vector<std::vector<float> > &wy, const std::vector<std::vector<float> > &wxy,
		const std::vector<std::vector<float> > &W3iji,
		const bool gcpflagin, const std::vector<int> &gcpinfoin, const std::vector<int> &gcptypein, const std::vector<std::vector<double> > &gcpdatain,
		const float ihw2i,
		const bool robustflagi,
		const bool itbreakflagi) : xg(xgin), ihw2(ihw2i), defdiff(5000. / sqrtf(ihw2i)) /* no more than 5 pixels */, wxi(wx), wyi(wy), wxyi(wxy)
{

	W3ij = W3iji;

	XXnr = XXnrin;

	int nrimagesin = (int) PMatricesi.size();
	PMatrices.resize(nrimagesin);
	for (int i = 0; i < nrimagesin; ++i)
		PMatrices[i] = PMatricesi[i];

	imgindflag = imgindflagin;

	if (imgindflag) {
		imgind = imgindin;
		nrimages = (int) imgind.size();
		}
	else {
		nrimages = nrimagesin;
		imgind.resize(nrimages);
		for (int i = 0; i < nrimages; ++i)
			imgind[i] = i;
	}

	minimagenr = 2;
	minimagenrdel = 2;
	if (nrimages > 2)
		minimagenrdel = 3;

	gcpflag = gcpflagin;
	aoflag = gcpflagin;

	PointXXvIn(XXi, pointXXvio, gcpinfoin);

	nrobsgcp = nrgcp = 0;

	// gcpinfo is fixed in pointXXvin!!
	if(gcpflag){
		for (int i = 0; i < (int) pointXXv.size(); ++i) {
			if(gcpinfo[i] > -1) {
				nrobsgcp += 3;
				++nrgcp;
			} // if gcpinfo
		} //
	} // if gcpflag

	if (aoflag) {
		nrobsao = nrobsgcp;
		nrao = nrgcp;
		if (robustflagi) {
			sig0aorv.resize(nrao);
			sig0aorvs.resize(nrao);
			W3ijold = W3ij;
			for (int i = 0; i < nrao; ++i) {
				sig0aorv[i] = sig0aorvs[i] = 0.f;
			} // for i
		}
	} // if aoflag

	PrepareCameras();

	PrepareIteration(cameras, addflagsin, gcptypein, gcpdatain);

	if (robustflagi) {
		sig00v.resize(nrimages);
		sig0rmvvs.resize(nrimages);
	}

	// XXnr is taken to be on the save side
	int xpmsize = nrimages * 6 + 3 * XXnr + 7 * nrcameras;
	XPMatrices.resize(xpmsize);
	XPMatricesnew.resize(xpmsize);

	return;
}



void LSA::PrepareIteration(const std::vector<Camera> &cset,
		const std::vector<bool> &addflagsin,
		const std::vector<int> &gcptypein, const std::vector<std::vector<double> > &gcpdatain) {

	addflags = addflagsin;

	nraddpi = 0;
	for (std::size_t i = 0; i < addflags.size(); ++i)
		if (addflags[i])
			++nraddpi;

	camparv.clear();
	camparv.resize(nrimages + 2);
	std::vector<int>::iterator pxi;
	for (pxi = camparv.begin(); pxi != camparv.end(); ++pxi)
		(*pxi) = 0;

	esvm.resize(3);

	eRM.resize(nrimages);
	ePC.resize(3,nrimages);
	eRMnew.resize(nrimages);
	ePCnew.resize(3,nrimages);

	const int imgind0 = imgind[0];
	RectifyRotationMatrix4(PMatrices[imgind0]);

	eRMnew[0] = PMatrices[imgind0].leftCols(3);
	ePCnew.col(0) = eRMnew[0].transpose() * PMatrices[imgind0].col(3);

	ex3u(2) = 1.;

	nrparir = 5 + nraddpi;
	nrpari = 6 + nraddpi;

	eBpi23.resize(2,nraddpi);
	eBpit23.resize(nraddpi,2);
	eBpi2r.resize(2,nrparir);
	eBpi2.resize(2,nrpari);
	eBpi2rt.resize(nrparir,2);
	eBpi2t.resize(nrpari,2);
	eBpi3r.resize(3,nrparir);
	eBpi3.resize(3,nrpari);
	eBpi3rt.resize(nrparir,3);
	eBpi3t.resize(nrpari,3);

	if(gcpflag) {
		if (gcpinfo.size() != pointXXv.size()) {
			std::cout << "Error with GCP: " << gcpinfo.size() << " " << pointXXv.size() << "\n";
			return;
		}
		const int gcpd = (int) gcpdatain.size();
		gcptype.resize(gcpd);
		gcpdata.resize(gcpd);
		for(int k = 0; k < gcpd; ++k) {
			gcptype[k] = gcptypein[k];
			gcpdata[k].resize(3);
			for(int j = 0; j < 3; ++j)
				gcpdata[k][j] = gcpdatain[k][j];
		} // for k
	} //if gcpflag

	nrcameras = 0;

	std::vector<int> cameracount(cset.size(), 0);

	nrcameras = (int) cset.size();
	cameras.resize(nrcameras);
	camerasnew.resize(nrcameras);
	camerasold.resize(nrcameras);
	for(int i = 0; i < nrcameras; ++i) {
		const Camera cam = cset[i];

		LSACamera& lcam = cameras[i];
		lcam.distpara.resize(2);
		lcam.K = cam.K.cast <double> ();
		lcam.distpara[0] = cam.distpara[0];
		lcam.distpara[1] = cam.distpara[1];

		camerasnew[i].K = lcam.K;
		camerasold[i].K = lcam.K;
		camerasnew[i].distpara = lcam.distpara;
		camerasold[i].distpara = lcam.distpara;
	}

	nrcameras = (int) cameras.size();


	lflflag = false;
	for(int i = 0; i < nrimages; ++i) {
		const Camera cam = cset[i];
		if (cam.K(0,0) > 2.f || cam.K(1,1) > 2.f) {
			lflflag = true;
			break;
			} // if
	} // for

	if (!aoflag)
		nrioff = nraddpi + nrparir + nrpari * (nrimages - 2);
	else
		nrioff = nrpari * nrimages;


	return;
} // LSA::PrepareIteration



void LSA::InputLSAP(const std::vector<bool> &addflagsin,
		const std::vector<Camera> &cset,
		const std::vector<int> &imgindin, std::vector<std::vector<IMAPNR> > &pointXXvio,
		const std::vector<PMatd> &PMatricesi,
		const Eigen::MatrixXd &XXi, const int XXnri,
		const std::vector<int> &gcpinfoin) {

	imgind = imgindin;
	nrimages = (int) imgind.size();

	int imgind0 = imgind[0];
	RectifyRotationMatrix4(PMatrices[imgind0]);

	eRMnew[0] = PMatrices[imgind0].leftCols(3);
	ePCnew.col(0) = eRMnew[0].transpose() * PMatrices[imgind0].col(3);

	gcpflag = aoflag = true;

	XXnr = XXnri;

	for (int i = 0; i < nrimages; ++i)
		PMatrices[imgind[i]] = PMatricesi[imgind[i]];

	PointXXvIn(XXi, pointXXvio, gcpinfoin);

	nrobsgcp = nrgcp = 0;

	// gcpinfo is fixed in pointXXvin!!
	for (int i = 0; i < (int) pointXXv.size(); ++i) {
		if (gcpinfo[i] > -1) {
			nrobsgcp += 3;
			++nrgcp;
		} // if gcpinfo
	} //

	nrobsao = nrobsgcp;
	nrao = nrgcp;

	if (gcpinfo.size() != pointXXv.size()) {
		std::cout << "Error with GCP: " << gcpinfo.size() << " " << pointXXv.size() << "\n";
		return;
	}

	addflags = addflagsin;

	nraddpi = 0;
	for (std::size_t i = 0; i < addflags.size(); ++i)
		if (addflags[i])
			++nraddpi;

	nrparir = 5 + nraddpi;
	nrpari = 6 + nraddpi;

	eBpi23.resize(2,nraddpi);
	eBpit23.resize(nraddpi,2);
	eBpi2r.resize(2,nrparir);
	eBpi2.resize(2,nrpari);
	eBpi2rt.resize(nrparir,2);
	eBpi2t.resize(nrpari,2);
	eBpi3r.resize(3,nrparir);
	eBpi3.resize(3,nrpari);
	eBpi3rt.resize(nrparir,3);
	eBpi3t.resize(nrpari,3);

	nrcameras = 0;

	std::vector<int> cameracount(cset.size(), 0);

	nrcameras = (int) cset.size();
	cameras.resize(nrcameras);
	camerasnew.resize(nrcameras);
	camerasold.resize(nrcameras);
	for(int i = 0; i < nrcameras; ++i) {
		const Camera cam = cset[i];

		LSACamera& lcam = cameras[i];
		lcam.distpara.resize(2);
		lcam.K = cam.K.cast <double> ();
		lcam.distpara[0] = cam.distpara[0];
		lcam.distpara[1] = cam.distpara[1];

		camerasnew[i].K = lcam.K;
		camerasold[i].K = lcam.K;
		camerasnew[i].distpara = lcam.distpara;
		camerasold[i].distpara = lcam.distpara;
	}

	nrcameras = (int) cameras.size();

	nrioff = nrpari * nrimages;

	int xpmsize = nrimages * 6 + 3 * XXnr + 7 * nrcameras;
	XPMatrices.resize(xpmsize);
	XPMatricesnew.resize(xpmsize);

	return;
} // InputLSAP



void LSA::Adjust(const bool robustflagi,
		std::ostream& os, const bool itbreakflagi, const int itmaxnumin, const std::vector<std::vector<float> > &W3ijit,
		const bool gcpflagin, const std::vector<int> &gcptypeit, const std::vector<std::vector<double> > &gcpdatait,
		const std::vector<Camera> &cset) {

	int itt, iter, iflag;
	bool ittflag;
	float sig0aoin = 1.e20f;

	sig0aoiold = sig0aoi = 1.e20f;

	int itmaxnum = 200;
	if (itmaxnumin > -1)
		itmaxnum = itmaxnumin;

	itt = 0;
	iniflag = true;
	ittflag = true;

	gcpflag = gcpflagin;
	aoflag = aoflag || gcpflagin;

	// sig0frac controls together with iflag when to stop the iteration.
	// If the fraction of old and new sig0 is smaller than sig0frac,
	// then there is not enough progress any more. As sometimes
	// there is an improvement only after a few iterations iflag > 0
	// avoids directly stopping if there is not an improvement for
	// an iteration

	float sig0frac = 1.05f; // 1.05f
	if ((!itbreakflagi && aoflag) || lflflag)
		sig0frac = 1.01f;

	double multfact = 1.e-3; // e-4

	if (itbreakflagi)
		multfact = 1.e-4;

	for (int it = 0; it < 2; ++it){

		sig0old = 1.e20f;

		if (it > 0) {
			iniflag = true;
			multfact = 1.e-6;
		}

		imagepointnr = 0;
		std::vector<std::vector<IMPNR> >::iterator pxv;
		for(pxv = pointXXv.begin(); pxv != pointXXv.end(); ++pxv)
			imagepointnr += (int) pxv->size();
		nrobs = 2 * imagepointnr;

		nrobsgcp = 0;
		nrgcp = 0;

		if(gcpflag){

			for(std::size_t i = 0; i < pointXXv.size(); ++i)
				if(gcpinfo[i] > -1) {
					nrobsgcp += 3;
					++nrgcp;
				}
		} //if

		nrobsao = nrobsgcp;
		nrobs += nrobsao;
		nrao = nrgcp;

		nrobv.resize(nrimages + 1);
		std::vector<int>::iterator pxi;
		for(pxi = nrobv.begin(); pxi != nrobv.end(); ++pxi)
			(*pxi) = 0;
		maxnrobsc = 0;
		for(pxv = pointXXv.begin(); pxv != pointXXv.end(); ++pxv) {
			if ((int) pxv->size() > maxnrobsc)
				maxnrobsc = (int) pxv->size();
			std::vector<IMPNR>::iterator pxvi;
			for(pxvi = pxv->begin(); pxvi != pxv->end(); ++pxvi)
				++nrobv[pxvi->iimage + 1];
		}

		for(int i = 1; i < nrimages+1; ++i) {
			if (nrobv[i] == 0){
				//----------------------------------------------------------------------------
				if (nrimages > 3)
				{
					std::cout << "0000000000000000000 " << i << '\n';
					os << "BA-Error: No observations in image ";
					os << i-1 << " (image index)\n";
				}
				return;
			} //if nrobv[i]
		}

		for(int i = 0; i < nrimages; ++i)
			nrobv[i+1] += nrobv[i];

		nrobvo = nrobv;

		if (ittflag) {
			Wij.resize(imagepointnr);
			ISflag.resize(imagepointnr);
			ISflagnew.resize(imagepointnr);
			for (int i = 0; i < imagepointnr; ++i)
				ISflag[i] = ISflagnew[i] = false;
			PrepareWeights2D(it, robustflagi);

			ittflag = false;

		} // if ittflag

		if (aoflag) {
			PrepareWeights3D(robustflagi,itbreakflagi,W3ijit);

			wsumg = (float) (2 * imagepointnr + nrobsao) / (wsumao + wsum);
			wsumao = (float) nrobsao / wsumao;
		}

		wsum = (float) (2 * imagepointnr) / wsum;

		if (!aoflag)
			wsumg = wsum;

		for(int j = 0; j < 3; ++j)
			for(int i = 0; i < XX.cols(); ++i)
				XPMatrices[nrioff + i * 3 + j] = XX(j,i);

		EncodeCameras();

		nrunknowns = nrioff3 + nraddpi * nrcameras;

		XXnew = XX;

		x.resize(nrunknowns);

		redundancy = nrobs - nrunknowns;

		if (redundancy < 0){
			if (nrimages > 3) {
				std::cout << "#######  Redundancy < 1  #########\n";
				os << "BA-Error: Low redundancy\n";
			}
			return;
		} // if redundancy

		// calc sig0old
		if (it == 0) {
			sig0rmv.resize(imagepointnr);
			DecodeCameras(XPMatrices);
			nrobv = nrobvo;
			sig0ao = 0.f;
			CalcCorrections(false, 1. + multfact, XPMatrices, true, false, robustflagi, itbreakflagi, W3ijit, gcptypeit, gcpdatait);

			if (aoflag) {
				sig0aoi = sig0 * wsumg / (float) redundancy;
				sig0aoin = sig0aoi;
				sig0 = (sig0 + sig0ao) * wsumg / (float) (redundancy + nrobsao - 7);
			}
			else {
				sig0 = sig0 * wsumg / (float) redundancy;
			}
			sig0in = sig0old = sig0;
		}
		else {
			sig0in = sig0old;
			sig0aoin = sig0aoiold;
		}

		// iflag is used to control robust adjustment
		iflag = 2;

		if (!robustflagi)
			iflag = 1;

		if (aoflag)
			iflag = 4;

		iter = 1;

		const int nrunknownsgiven = nrunknowns;
		const int nrobsgiven = nrobs;

		while(iflag > 0){

			if (it == 0 && robustflagi) {

				nrunknowns = nrunknownsgiven;
				nrobs = nrobsgiven;

				Reweight2D(iflag);

				redundancy = nrobs - nrunknowns;

				wsum = (float) (2 * imagepointnr) / wsum;

				if (!aoflag)
					wsumg = wsum;
			}

			if (CalcBundleSep(multfact, robustflagi, itbreakflagi, W3ijit, gcptypeit, gcpdatait)) {
				if (nrimages > 3)
					os << "BA-Error: CalcBundleSep\n";
				return;
			}


			//improved value for vector of unknowns

			for(int i = 0; i < nrunknowns; ++i) {
				if (boost::math::isnan(x[i]) || boost::math::isinf(x[i])){
					std::cout << "NANNANNANXPM\n";
					std::cout << "BA-Error: NANs in result vector " << i << " " << x[i] << " " << nrioff << " " << nrioff3 << "\n";
					return;
				} //if isnan x
				XPMatricesnew[i] = XPMatrices[i] - x[i];
			}

			DecodeCameras(XPMatricesnew);
			nrobv = nrobvo;
			CalcCorrections(false, 1. + multfact, XPMatricesnew, true, false, robustflagi, itbreakflagi, W3ijit, gcptypeit, gcpdatait);

			if (aoflag) {
				sig0aoi = sig0 * wsumg / (float) redundancy;
				sig0 = (sig0 + sig0ao) * wsumg / (float) (redundancy + nrobsao - 7);
			}
			else
				sig0 = sig0 * wsumg / (float) redundancy;


			if (sig0old < sig0 * sig0frac)
				--iflag;

			if (sig0 <= sig0old || iter == 1) { //  || itt == 0) {
				StoreResultsIteration(itbreakflagi);

				sig0old = sig0;
				sig0aoiold = sig0aoi;
			}

			if (itbreakflagi && iter > itmaxnum) // 0
				break;

			if (itt > 200) // 200
				break;

			++itt;
			++iter;

		} //while iflag > 0

		const float _sig0 = sqrtf(sig0in * ihw2);
		const float _sig1 = sqrtf(sig0old * ihw2);

		os << "Sigma0: ";
		if (aoflag)
			os << sqrtf(sig0aoin * ihw2) << ' ' << sqrtf(sig0aoiold * ihw2) << "  ";
		else
			os << _sig0 << ' ' << _sig1 << "  ";
		if (it == 0 && robustflagi) {
			os << redundancynew << ' ' << nrobsnew << ' ' << nrunknownsnew;
		}
		else {
			os << redundancy << ' ' << nrobs << ' ' << nrunknowns;
		}
		os << "  " << "IT: " << it << "  ITT: " << itt << '\n';
		if (aoflag) {
			os << "Sigma0AO: " << sqrtf(sig0ao * wsumao / (float) ((nrobsao - 7))) << "  " << _sig0 << ' ' << _sig1 << "  " << nrobsao << " " << wsumao << '\n';
		}

		if (itbreakflagi)
			return;

		// Throw out outliers; criterion: weight in x or y direction > 4
		if (robustflagi && it == 0) {
			RobustDelete();
		}

		nrobv = nrobvo;
		ittflag = true;

	} // for it

	return;
} //LSA::Adjust



void LSA::OutputLSAP(std::vector<PMatd> &PMatricesi,
		std::vector<std::vector<IMAPNR> > &pointXXvio,Eigen:: MatrixXd &XXio,
		std::vector<int> &gcpinfoo, std::vector<Camera> &cameraso) {

	int nrimagesin = (int) PMatricesi.size();

	PointXXvOut(XXio, pointXXvio, gcpinfoo);

	//restore values of optimum!!!!
	for (int i = 0; i < nrimagesin; ++i)
		PMatricesi[i] = PMatrices[i];

	for(int i = 0; i < nrcameras; ++i){
		cameraso[i] = Camera();
		cameraso[i].K = camerasnew[i].K.cast <float> ();
		cameraso[i].distpara[0] = (float) camerasnew[i].distpara[0];
		cameraso[i].distpara[1] = (float) camerasnew[i].distpara[1];
	} //for k

	return;
} // OutputLSAP



bool LSA::OutputLSAc(std::vector<PMat> &PMatricesi,
		std::vector<std::vector<IMAPNR> > &pointXXvio, Eigen::MatrixXf &XXio,
		std::vector<Camera> &cameraso) {

	std::vector<int> gcpinfoo;

	Eigen::MatrixXd XXX = XXio.cast <double> ();
	PointXXvOut(XXX, pointXXvio,gcpinfoo);

	XXio = XXX.cast <float> ();

	for (int i = 0; i < nrimages; ++i)
		PMatricesi[i] = PMatrices[i].cast <float> ();

	for(int i = 0; i < nrcameras; ++i){
		cameraso[i] = Camera();
		cameraso[i].K = camerasnew[i].K.cast <float> ();
		cameraso[i].distpara[0] = (float) camerasnew[i].distpara[0];
		cameraso[i].distpara[1] = (float) camerasnew[i].distpara[1];
	} //for i

	return true;
} // OutputLSAc



void LSA::PrepareCameras() {

	Eigen::MatrixXd Tv(3,nrimages);

	for (int i = 0; i < nrimages; ++i){
		const int imgindi = imgind[i];
		Tv.col(i) = PMatrices[imgindi].leftCols(3) * PMatrices[imgindi].col(3);
	}

	// T0 is vector to first camera
	Eigen::Vector3d T0 = Tv.col(0);


	double normdist = 0.f;
	double dist;
	for (int i = 1; i < nrimages; ++i){
		dist = 0;
		for (int j = 0; j < 3; ++j)
			dist += (Tv(j,i) - T0(j)) * (Tv(j,i) - T0(j));
		if (dist > normdist){
			normdist = dist;
			normcamera = i;
		}
	}

	// To avoid that one image is set to "normcamera"
	if (aoflag)
		normcamera = -1;

	return;
}



void LSA::PrepareWeights2D(const int it, const bool robustflagi) {

	if (robustflagi) {
		wsum = 0.f;

		int i;

		if (it == 0) {
			Wijoldv.resize((int) pointXXv.size());
			std::vector<std::vector<IMPNR> >::iterator pxv;
			std::vector<std::vector<std::vector<float> > >::iterator wxyo;
			for(pxv = pointXXv.begin(), i = 0, wxyo = Wijoldv.begin(); pxv != pointXXv.end(); ++pxv, ++i, ++wxyo){
				wxyo->resize((int) pxv->size());
				std::vector<IMPNR>::iterator pxvi;
				std::vector<std::vector<float> >::iterator wxyoi;
				for(pxvi = pxv->begin(), wxyoi = wxyo->begin(); pxvi != pxv->end(); ++pxvi, ++wxyoi){
					const int jjj = imgind[pxvi->iimage];
					wxyoi->resize(3);
					const int pxx = pxvi->pointnr;
					(*wxyoi)[0] = wxi[jjj][pxx];
					(*wxyoi)[1] = wyi[jjj][pxx];
					(*wxyoi)[2] = wxyi[jjj][pxx];

					wsum += (*wxyoi)[0] + (*wxyoi)[1];
				} // for pxvi
			} //for pxv
		} // if it
		else {
			std::vector<std::vector<IMPNR> >::iterator pxv;
			std::vector<std::vector<std::vector<float> > >::iterator wxyo;
			for(wxyo = Wijoldv.begin(), pxv = pointXXv.begin(); wxyo != Wijoldv.end(); ++wxyo, ++pxv) {
				std::vector<std::vector<float> >::iterator wxyoi;
				for(wxyoi = wxyo->begin(); wxyoi != wxyo->end(); ++wxyoi)
					wsum += ((*wxyoi)[0] + (*wxyoi)[1]);
			} // for wxy0
		}

		nrobv = nrobvo;

		std::vector<std::vector<IMPNR> >::iterator pxv;
		std::vector<std::vector<std::vector<float> > >::iterator wxyo;
		for(pxv = pointXXv.begin(), wxyo = Wijoldv.begin(), i = 0; pxv != pointXXv.end(); ++pxv, ++wxyo, ++i) {
			std::vector<IMPNR>::iterator pxvi;
			std::vector<std::vector<float> >::iterator wxyoi;
			for(pxvi = pxv->begin(), wxyoi = wxyo->begin(); pxvi != pxv->end(); ++pxvi, ++wxyoi) {
				const int jjj = pxvi->iimage;
				const int nrjjj = nrobv[jjj];
				Wij[nrjjj].resize(3);
				Wij[nrjjj][0] = (*wxyoi)[0];
				Wij[nrjjj][1] = (*wxyoi)[1];
				Wij[nrjjj][2] = (*wxyoi)[2];
				++nrobv[jjj];
			} //for pxvi
		} // for pxv
	} // if (robustflagi
	else {
		nrobv = nrobvo;

		wsum = 0.f;

		std::vector<std::vector<IMPNR> >::iterator pxv;
		for(pxv = pointXXv.begin(); pxv != pointXXv.end(); ++pxv) {
			std::vector<IMPNR>::iterator pxvi;
			for(pxvi = pxv->begin(); pxvi != pxv->end(); ++pxvi) {
				const int jjj = pxvi->iimage;
				const int pxx = pxvi->pointnr;
				const int nrjjj = nrobv[jjj];
				const int ijjj = imgind[jjj];
				Wij[nrjjj].resize(3);
				Wij[nrjjj][0] = wxi[ijjj][pxx];
				Wij[nrjjj][1] = wyi[ijjj][pxx];
				Wij[nrjjj][2] = wxyi[ijjj][pxx];
				wsum += wxi[ijjj][pxx] + wyi[ijjj][pxx];
				++nrobv[jjj];
			} //for pxvi
		} // for pxv
	} // else

	return;
} // PrepareWeights2D



void LSA::PrepareWeights3D(const bool robustflagi,
		const bool itbreakflagi, const std::vector<std::vector<float> > &W3ijit) {

	nrobsgcp = nrgcp = 0;

	if (gcpflag)
		for (int i = 0; i < (int) gcpinfo.size(); ++i)
			if (gcpinfo[i] > -1) {
				nrobsgcp += 3;
				++nrgcp;
			}

	nrobsao = nrobsgcp;
	nrao = nrgcp;

	if (robustflagi)
		W3ij = W3ijold;

	wsumao = 0.f;

	if (itbreakflagi) {
		for (int i = 0; i < nrao; ++i)
			wsumao += W3ijit[0][i] + W3ijit[1][i] + W3ijit[2][i];
	}
	else {
		for (int i = 0; i < nrao; ++i) {
			wsumao += W3ij[0][i] + W3ij[1][i] + W3ij[2][i];
		}
	}

	return;
} // PrepareWeights3D



// This implements re-weighting.
void LSA::Reweight2D(const int iflag) {

	for (int i = 0; i < nrimages; ++i)
		sig0rmvvs[i].resize(0);

	nrobv = nrobvo;
	std::vector<std::vector<IMPNR> >::iterator pxv;
	for(pxv = pointXXv.begin(); pxv != pointXXv.end(); ++pxv) {
		int pxs = (int) pxv->size() * 2;
		if (pxs < 4)
			pxs = 4;
		const float rf = (float) pxs / (float) (pxs - 3);
		const float redfac = rf * rf;
		std::vector<IMPNR>::iterator pxvi;
		for(pxvi = pxv->begin(); pxvi != pxv->end(); ++pxvi){
			const int kkk = pxvi->iimage;
			sig0rmvvs[kkk].push_back(sig0rmv[nrobv[kkk]] * redfac);
			++nrobv[kkk];
		} // for pxvi
	} // for pxv
	for (int i = 0; i < nrimages; ++i) {
		std::sort(sig0rmvvs[i].begin(), sig0rmvvs[i].end());
		sig00v[i] = 1.f / (sig0rmvvs[i][(int) sig0rmvvs[i].size() / 2] * 2.1981f); // 1.4826 // 1.706f - 1.306^2
		if (iflag > 1)
			sig00v[i] *= 0.001f;
	}

	nrobv = nrobvo;
	wsum = 0.f;

	std::vector<std::vector<std::vector<float> > >::iterator wxyo;
	for(pxv = pointXXv.begin(), wxyo = Wijoldv.begin(); pxv != pointXXv.end(); ++pxv, ++wxyo){
		std::vector<IMPNR>::iterator pxvi;
		std::vector<std::vector<float> >::iterator wxyoi;
		for(pxvi = pxv->begin(), wxyoi = wxyo->begin(); pxvi != pxv->end(); ++pxvi, ++wxyoi){
			const int kkk = pxvi->iimage;
			const int nro = nrobv[kkk];

			float ww = sig0rmv[nro] * sig00v[kkk];

			ISflag[nro] = ww > kf;

			if (ISflag[nro]) {
				ww = 1.e-4f;
				nrobs -= 2;
			}
			else {
				ww = 1.f;
			}

			Wij[nro][0] = ww * (*wxyoi)[0];
			Wij[nro][1] = ww * (*wxyoi)[1];
			Wij[nro][2] = ww * (*wxyoi)[2];

			if (Wij[nro][0] + Wij[nro][1] < 2.e-8f) {
				Wij[nro][0] = Wij[nro][1] = 1.e-8f;
				Wij[nro][2] = 0.f;
			}
			wsum += Wij[nro][0] + Wij[nro][1];
			++nrobv[kkk];
		} //for pxvi
	} //for pxv

	return;
} // Reweight2D



void LSA::StoreResultsIteration(const bool itbreakflagi) {

	XPMatrices = XPMatricesnew;

	for(int k = 0; k < nrimages; ++k) {
		const int imgk = imgind[k];

		PMatrices[imgk].leftCols(3) = eRMnew[k];
		PMatrices[imgk].col(3) = eRMnew[k] * ePCnew.col(k);
	}


	for (int k = 0; k < nrcameras; ++k) {
		cameras[k].K = camerasnew[k].K;
		cameras[k].distpara = camerasnew[k].distpara;
	}

	XX = XXnew;


	ISflagnew = ISflag;

	nrobsnew = nrobs;
	nrunknownsnew = nrunknowns;
	redundancynew = redundancy;

	return;
} // StoreResultsIteration



void LSA::RobustDelete() {

	int i, k, ks, psdiff;
	bool pflag;

	nrobv = nrobvo;
	erasevpts.resize(0);
	std::vector<std::vector<IMPNR> >::iterator pxv;
	for(pxv = pointXXv.begin(), i = 0; pxv != pointXXv.end(); ++pxv, ++i){
		pflag = false;
		erasev.resize(0);
		std::vector<IMPNR>::iterator pxvi;
		for(pxvi = pxv->begin(), ks = 0; pxvi != pxv->end(); ++pxvi, ++ks){
			k = pxvi->iimage;
			if (ISflagnew[nrobv[k]]) {
				erasev.push_back(ks);
				pflag = true;
			} //if ISflagnew
			++nrobv[k];
		} // for pxvi
		if (pflag){
			psdiff = (int) pxv->size() - (int) erasev.size();
			if ((!gcpflag && psdiff < minimagenrdel) || (gcpflag && gcpinfo[i] > -1 && psdiff < 1) ||
					(gcpflag && gcpinfo[i] < 0 && psdiff < minimagenrdel)){
				pointXXv[i].clear();
				erasevpts.push_back(i);
			} // if
			else{
				for (k = (int) erasev.size()-1; k > -1; --k){
					pxv->erase(pxv->begin()+erasev[k]);
					Wijoldv[i].erase(Wijoldv[i].begin()+erasev[k]);
				} // for k
			} // else
		} // if pflag
	} // for pxv


	sort(erasevpts.begin(),erasevpts.end());

	const bool gcpf = gcpflag || (int) gcpinfo.size() > 0;

	for(i = (int) erasevpts.size() - 1; i > -1; --i){
		const int evs = erasevpts[i];
		pointXXv.erase(pointXXv.begin() + evs);
		nrpv.erase(nrpv.begin() + evs);
		Wijoldv.erase(Wijoldv.begin() + evs);
		if (gcpf)
			gcpinfo.erase(gcpinfo.begin() + evs);
	}
	nrpts -= (int) erasevpts.size();

	if ((int) erasevpts.size() > 0) {
		Eigen::MatrixXd XXX(3, (int) XX.cols() - (int) erasevpts.size());

		k = 0;
		int j = 0;
		for (i = 0; i < XX.cols(); ++i) {
			if (i == erasevpts[k]) {
				if (k < (int) erasevpts.size() - 1)
					++k;
				continue;
			}
			XXX(0,j) = XX(0,i);
			XXX(1,j) = XX(1,i);
			XXX(2,j) = XX(2,i);
			++j;
		} // for i

		XX = XXX;
		XXnew = Eigen::MatrixXd (3, XX.cols());
	}

	return;
}



void LSA::EncodeCameras() {
	int nrimagesall, nrimagesoffset;

	for(int k = 0; k < nrimages + 1; ++k)
		camparv[k] = 0;

	if (aoflag) {
		nrimagesall = nrimages;
		camparv[1] = nrpari;
		nrimagesoffset = 0;
	}
	else {
		camparv[1] = nraddpi;
		nrimagesall = nrimages - 1;
		nrimagesoffset = 1;
	}

	if (!aoflag) {
		int addp = -1;
		const int imgind0 = imgind[0];
		if (addflags[0])
			XPMatrices[++addp] = cameras[imgind0].distpara[0];
		if (addflags[1])
			XPMatrices[++addp] = cameras[imgind0].distpara[1];
		if (addflags[2])
			XPMatrices[++addp] = cameras[imgind0].K(0,0);
	}

	for (int k = 0; k < nrimagesall; ++k) {

		const int imgindkn = imgind[k + nrimagesoffset];
		const int campkn = camparv[k + nrimagesoffset];

		RectifyRotationMatrix4(PMatrices[imgindkn]);

		eRM[k] = PMatrices[imgindkn].leftCols(3);

		ePC.col(k) = eRM[k].transpose() * PMatrices[imgindkn].col(3);

		for(int i = 0; i < 3; ++i)
			XPMatrices[campkn + i] = 0.;


		if (k == normcamera - 1) { // normcamera (k == 0)
			float maxv = 0.f;
			int maxi = 0;
			for(int i = 0; i < 3; ++i) {
				if(fabs(ePC(i,k)) > maxv){
					maxv = (float) fabs(ePC(i,k));
					maxi = i;
				} //if abs
			} //for i

			int ii = -1;
			for(int i = 0; i < 3; ++i) {
				if(!(maxi == i))
					esvm[++ii] = i;
			} //for i
			esvmaxm = maxi;

			for(int i = 0; i < 2; ++i)
				XPMatrices[campkn + i + 3] = 0.;
			int addp = 4 + campkn;
			if (addflags[0])
				XPMatrices[++addp] = cameras[imgindkn].distpara[0];
			if (addflags[1])
				XPMatrices[++addp] = cameras[imgindkn].distpara[1];
			if (addflags[2])
				XPMatrices[++addp] = cameras[imgindkn].K(0,0);
			camparv[k + 2] = campkn + nrparir;
		} //if k ==
		else{
			for(int i = 0; i < 3; ++i)
				XPMatrices[campkn + i + 3] = 0.;
			int addp = 5 + campkn;
			if (addflags[0])
				XPMatrices[++addp] = cameras[imgindkn].distpara[0];
			if (addflags[1])
				XPMatrices[++addp] = cameras[imgindkn].distpara[1];
			if (addflags[2])
				XPMatrices[++addp] = cameras[imgindkn].K(0,0);
			camparv[k + 2] = camparv[k + 1] + nrpari;
		} //else
	} //for k nrimages-1

	nrioff3 = nrioff + 3 * nrpts;

	return;
}



// Bundle solution separated for 3D points, camera matrices and
// additional parameters. The algorithm basically follows the book
// of Mikhail et al. on "Introduction to Modern Photogrammetry".
bool LSA::CalcBundleSep(const double multfact, const bool robustflagi,
		const bool itbreakflagi, const std::vector<std::vector<float> > &W3ijit,
		const std::vector<int> &gcptypeit, const std::vector<std::vector<double> > &gcpdatait) {

	int nbegin, nend = 0;

	int nrimagesapc = nrimages;

	sig0ao = 0.f;

	begini = 1;
	if (aoflag)
		begini = 0;

	const double multfact1 = 1. + multfact;

	if (iniflag) {
		exiv.resize(nrimagesapc);
		Npqqnumv.resize(nrimagesapc);
		for (int i = 0; i < nrimagesapc; ++i) {
			Npqqnumv[i].resize(0);
			Npqqnumv[i].push_back(true);
		}
		eNpqqv.resize(nrimagesapc);
		for(int i = 0; i < nrimages; ++i)
			eNpqqv[i].resize(1);
		colv.resize(nrimages);
		for(int i = begini; i < nrimages; ++i)
			colv[i] = i;

		nbeginv.resize(nrimages + 2);
		nendv.resize(nrimages + 2);
	} // if iniflag


	// Preparation loop for images
	for(int i = begini; i < nrimages; ++i){
		nbeginv[i] = camparv[i];
		nendv[i] = camparv[i + 1];

		nbegin = nbeginv[i];
		nend = nendv[i];

		if (iniflag) {
			eNpqqv[i][0].resize(nend - nbegin, nend - nbegin);
			exiv[i].resize(nend - nbegin, 1);
		}

		eNpqqv[i][0].setZero();
		exiv[i].setZero();

	}


	if (iniflag) {
		eNlij.resize(maxnrobsc);
	} // if iniflag

	for (int i = begini; i < nrimages; ++i){
		for (int j = 1; j < (int) eNpqqv[i].size(); ++j) {
			Npqqnumv[i][j] = false;
			eNpqqv[i][j].setZero();
		} // for (j
	} // for i


	nrobv = nrobvo;
	CalcCorrections(false, multfact1, XPMatrices, false, true, robustflagi, itbreakflagi, W3ijit, gcptypeit, gcpdatait);


	for(int i = begini; i < nrimagesapc; ++i)
		for(int j = 0; j < eNpqqv[i][0].rows(); ++j)
			eNpqqv[i][0](j,j) *= multfact1;

	for(int i = begini; i < nrimagesapc; ++i)
		for(int j = 0; j < eNpqqv[i][0].rows(); ++j)
			eNpqqv[i][0](j,j) += 1.e-2;

	int beginiold = begini;
	if (aoflag)
		begini = 0;

	if (SparsePCG(begini))
		return (true);

	begini = beginiold;

	for(int j = begini; j < nrimages; ++j){
		nbegin = nbeginv[j];
		nend = nendv[j];
		for(int k = nbegin; k < nend; ++k)
			x[k] = (float) exiv[j](k-nbegin);
	} //for j

	nrobv = nrobvo;
	CalcCorrections(true, multfact1, XPMatrices, false, true, robustflagi, itbreakflagi, W3ijit, gcptypeit, gcpdatait);


	iniflag = false;

	return (false);
} //LSA::CalcBundleSep



// Computation of points and camera matrices from given parameters
// leading to backprojection errors in the image plane
void LSA::CalcCorrections(const bool solveflag, const double multfact1, const std::vector<double> &XM, const bool srmvflag,
		const bool derivflag, const bool robustflagi,
		const bool itbreakflagi, const std::vector<std::vector<float> > &W3ijit,
		const std::vector<int> &gcptypeit, const std::vector<std::vector<double> > &gcpdatait) {

	sig0 = 0.f;

	int kkk = nrioff - 1;
	int gcpnum = 0;


	int i;
	std::vector<std::vector<IMPNR> >::iterator pxv;
	std::vector<std::vector<float> >::iterator wxyoi;
	for(pxv = pointXXv.begin(), i = 0; pxv != pointXXv.end(); ++pxv, ++i) {

		if (derivflag) {
			eNppi.setZero();
			etppi.setZero();
		} // if derivflag

		const int nrioff3i = nrioff + 3 * i;

		XXnew(0,i) = XM[nrioff3i];
		XXnew(1,i) = XM[nrioff3i + 1];
		XXnew(2,i) = XM[nrioff3i + 2];

		if (robustflagi)
			wxyoi = Wijoldv[i].begin();
		int ks = 0;
		for(std::vector<IMPNR>::iterator pxvi = pxv->begin(); pxvi != pxv->end(); ++pxvi, ++ks){
			const int k = pxvi->iimage;

			const int nro = nrobv[k];

			const int gp = pxvi->pointnr;

			ex3 = eRMnew[k] * (XXnew.col(i) + ePCnew.col(k));

			double xs3 = ex3(2);

			if (fabs(xs3) < 1.e-3) {
				if (xs3 < 0.)
					xs3 = -1.e3;
				else
					xs3 = 1.e3;
			}
			else
				xs3 = 1. / xs3;

			ex3(0) *= xs3;
			ex3(1) *= xs3;
			ex3(2) = 1.;

			double diff1, diff2;

			const int iik = imgind[k];
			int camerak = 0;

			camerak = imgind[k];
			ex3u(0) = ex3(0);
			ex3u(1) = ex3(1);
			RadialReconstructionImage(ex3u(0), ex3u(1), camerasnew[camerak]);

			const Eigen::Matrix3d &K = camerasnew[camerak].K;
			ex3c(0) = ex3u(0) * K(0,0) + ex3u(1) * K(0,1) + K(0,2);
			ex3c(1) = ex3u(1) * K(1,1) + K(1,2);

			// actual difference between given points and predicted points
			diff1 = (double) xg[iik](0,gp) - ex3c(0);
			diff2 = (double) xg[iik](1,gp) - ex3c(1);

			if (diff1 > defdiff) {
				std::cout << "DDD111: " << k << "  " << diff1 << "   " << defdiff << " " << '\n';
				diff1 = defdiff;
			}
			if (diff1 < -defdiff) {
				std::cout << "DDD111: " << k << "  " << diff1 << "   " << defdiff << " " << '\n';
				diff1 = -defdiff;
			}

			if (diff2 > defdiff) {
				std::cout << "DDD222: " << k << "  " << diff2 << "   " << defdiff << " " << '\n';
				diff2 = defdiff;
			}
			if (diff2 < -defdiff) {
				std::cout << "DDD222: " << k << "  " << diff2 << "   " << defdiff << " " << '\n';
				diff2 = -defdiff;
			}

			bool nanflag = false;
			if (boost::math::isnan(diff1) || boost::math::isinf(diff1)){
				std::cout << "NANX: " << XXnew.col(i).transpose() << ' ' << eRMnew[k] << '\n';
				nanflag = true;
			}

			if (boost::math::isnan(diff2) || boost::math::isinf(diff2)) {
				std::cout << "NANY: " << XXnew.col(i).transpose() << ' ' << eRMnew[k] << '\n';
				nanflag = true;
			}

			if (nanflag) {
				diff1 = diff2 = 0.;
				if (derivflag)
					CalcDerivatives(solveflag, camerak, k, ks, nro, xs3, diff1, diff2);
				sig0rmv[nro] = 1.e20f;
			} // if ISflag
			else {
				if (derivflag)
					CalcDerivatives(solveflag, camerak, k, ks, nro, xs3, diff1, diff2);

				const float diff1f = (float) diff1;
				const float diff2f = (float) diff2;
				const float d11 = diff1f * diff1f;
				const float d12 = 2.f * diff1f * diff2f;
				const float d22 = diff2f * diff2f;
				float sig0val,sig0rval;
				sig0val = (d11 * Wij[nro][0] + d12 * Wij[nro][2] + d22 * Wij[nro][1]);
				if (robustflagi) {
					sig0rval = d11 * (*wxyoi)[0] + d12 * (*wxyoi)[2] + d22 * (*wxyoi)[1];
					++wxyoi;
				}
				else
					sig0rval = sig0val;
				if (srmvflag)
					sig0rmv[nro] = sig0rval;
				sig0 += sig0val;
			} // else ISflag

			++nrobv[k];
		} //for pxvi

		if (gcpflag && gcpinfo[i] > -1) {
			if (itbreakflagi)
				CalcDerivatives3DGCP(i, gcpnum, XM, robustflagi, itbreakflagi, W3ijit, gcptypeit, gcpdatait, derivflag);
			else
				CalcDerivatives3DGCP(i, gcpnum, XM, robustflagi, itbreakflagi, W3ij, gcptypeit, gcpdatait, derivflag);
			++gcpnum;
		}

		if (derivflag) {
			for(int j = 0; j < 3; ++j)
				eNppi(j,j) *= multfact1;
			eNppi1.noalias() = eNppi.inverse();
			if (eNppi1.determinant() < 0.)
				std::cout << "IIIIIIIIIIIIIIIIIIIIIIII: " << i << " " << eNppi1.determinant() << " " << eNppi.determinant() << " " << multfact1 - 1. << '\n';

		} // if derivflag

		if (!solveflag && derivflag) {

			int is = 0;
			for(std::vector<IMPNR>::iterator pxvi = pxv->begin(); pxvi != pxv->end(); ++pxvi, ++is) {
				const int ij = pxvi->iimage;
				if (ij < begini)
					continue;

				eNlijN = eNlij[is] * eNppi1;
				eNpqqv[ij][0].noalias() -= eNlijN * eNlij[is].transpose();
				exiv[ij].noalias() -= eNlijN * etppi;
				// loop for images (columns)
				int iis;
				std::vector<IMPNR>::iterator pxvii;
				for(pxvii = pxvi + 1, iis = is + 1; pxvii != pxv->end(); ++pxvii, ++iis){
					const int ii = pxvii->iimage; // column of matrix
					if (ii < begini)
						continue;
					if (ij < ii) {
						if (ij < colv[ii]){
							colv[ii] = ij;
							for (int jj = (int) eNpqqv[ii].size(); jj < (ii - ij + 1); ++jj){
								eNpqqv[ii].push_back(Eigen::MatrixXd::Zero(nendv[ii-jj] - nbeginv[ii-jj], nendv[ii] - nbeginv[ii]));
								Npqqnumv[ii].push_back(false);
							} // for jj
						} // if ij
						eNpqqv[ii][ii - ij].noalias() -= eNlijN * eNlij[iis].transpose();
						Npqqnumv[ii][ii - ij] = true;
					} // if (ij < ii)
					else {
						if (ii < colv[ij]){
							colv[ij] = ii;
							for (int jj = (int) eNpqqv[ij].size(); jj < (ij - ii + 1); ++jj){
								eNpqqv[ij].push_back(Eigen::MatrixXd::Zero(nendv[ij-jj] - nbeginv[ij-jj], nendv[ij] - nbeginv[ij]));
								Npqqnumv[ij].push_back(false);
							} // for jj
						} // if ii < colv
						eNpqqv[ij][ij - ii].noalias() -= eNlij[iis] * eNlijN.transpose();
						Npqqnumv[ij][ij - ii] = true;
					} // else ij < ii
				} //for pxvii
			} //for pxvi
		} // if derivflag


		if (solveflag) {
			Eigen::Vector3d expp;
				expp = eNppi1 * etppi;

				int is;
				std::vector<IMPNR>::iterator pxvi;
				for(pxvi = pxv->begin(), is = 0; pxvi != pxv->end(); ++pxvi, ++is){
					const int j = pxvi->iimage;
					if (j >= begini) {
						expp.noalias() -= eNppi1 * eNlij[is].transpose() * exiv[j];
					}
				}
				x[++kkk] = expp(0);
				x[++kkk] = expp(1);
				x[++kkk] = expp(2);
		} // if solveflag
	} // for i

	return;
} //void LSA::CalcCorrections




// Computation of derivatives
void LSA::CalcDerivatives(const bool solveflag, const int camerak,
		const int k, const int ks, const int nro, const double xs3,
		const double diff1, const double diff2) {

	const double x1 = ex3(0) * xs3;
	const double x2 = ex3(1) * xs3;

	Eigen::Matrix<double, 2 , 3> eBpp;
	Eigen::Matrix<double, 3 , 2> eBppt;
	Eigen::Matrix<double, 2 , 5> eBpi5;
	Eigen::Matrix<double, 5 , 2> eBpi5t;
	Eigen::Matrix<double, 2 , 6> eBpi6;
	Eigen::Matrix<double, 6 , 2> eBpi6t;
	Eigen::Vector2d efp2(diff1, diff2);

	if (k >= begini) {
		const double xs31 = 1. / xs3;
		const double xs32 = xs31 * xs31;

		if (k == normcamera) {
			for (int j = 0; j < 2; ++j) {
				const int esvmj = esvm[j];
				eBpi2r(0, j + 3) = -eRMnew[k](0,esvmj) * xs3 + eRMnew[k](2,esvmj) * x1;
				eBpi2r(1, j + 3) = -eRMnew[k](1,esvmj) * xs3 + eRMnew[k](2,esvmj) * x2;
			}
			eBpi2r(0,0) = x2 * x1 * xs32;
			eBpi2r(1,0) = 1. + x2 * x2 * xs32;
			eBpi2r(0,1) = -1. - x1 * x1 * xs32;
			eBpi2r(1,1) = -x2 * x1 * xs32;
			eBpi2r(0,2) = x2 * xs31;
			eBpi2r(1,2) = -x1 * xs31;

			eCalcCalibratedDerivativesX(eBpi2r, camerak);
		} // if k == normcamera
		else {
			for(int j = 0; j < 3; ++j) {
				eBpi2(0, j + 3) = -eRMnew[k](0,j) * xs3 + eRMnew[k](2,j) * x1;
				eBpi2(1, j + 3) = -eRMnew[k](1,j) * xs3 + eRMnew[k](2,j) * x2;
			} // for j
			eBpi2(0,0) = x2 * x1 * xs32;
			eBpi2(1,0) = 1. + x2 * x2 * xs32;
			eBpi2(0,1) = -1. - x1 * x1 * xs32;
			eBpi2(1,1) = -x2 * x1 * xs32;
			eBpi2(0,2) = x2 * xs31;
			eBpi2(1,2) = -x1 * xs31;

			eCalcCalibratedDerivativesX(eBpi2, camerak);
		}

		if (k == normcamera) {
			for(int j = 0; j < 3; ++j){
				eBpp(0,j) = -eRMnew[k](0,j) * xs3 + eRMnew[k](2,j) * x1;
				eBpp(1,j) = -eRMnew[k](1,j) * xs3 + eRMnew[k](2,j) * x2;
			} //for j

			eCalcCalibratedDerivatives3(eBpp, camerak);
		}
		else {
			for(int j = 0; j < 3; ++j){
				eBpp(0,j) = eBpi2(0,j + 3);
				eBpp(1,j) = eBpi2(1,j + 3);
			}
		}
	} // if k >= begini


	if (k == 0 && (begini != 0 || (!aoflag))) {
		for(int j = 0; j < 3; ++j){
			eBpp(0,j) = -eRMnew[k](0,j) * xs3 + eRMnew[k](2,j) * x1;
			eBpp(1,j) = -eRMnew[k](1,j) * xs3 + eRMnew[k](2,j) * x2;
		} //for j

			eCalcCalibratedDerivatives3(eBpp, camerak);
	}

	const double xp2 = ex3(0) * ex3(0);
	const double yp2 = ex3(1) * ex3(1);
	const double xyp2 = xp2 + yp2;
	const double xyp4 = xyp2 * xyp2;
	if ((aoflag) || k > 0) {
		if (k == normcamera) {
			int addp = 4;
			if (addflags[0]) {
				++addp;
				eBpi2r(0,addp) = -ex3(0) * xyp2 * camerasnew[camerak].K(0,0) - ex3(1) * xyp2 * camerasnew[camerak].K(0,1);
				eBpi2r(1,addp) = -ex3(1) * xyp2 * camerasnew[camerak].K(1,1);
			}
			if (addflags[1]) {
				++addp;
				eBpi2r(0,addp) = -ex3(0) * xyp4 * camerasnew[camerak].K(0,0) - ex3(1) * xyp4 * camerasnew[camerak].K(0,1);
				eBpi2r(1,addp) = -ex3(1) * xyp4 * camerasnew[camerak].K(1,1);
			}
			if (addflags[2]) {
				++addp;
				eBpi2r(0,addp) = -ex3u(0);
				eBpi2r(1,addp) = -ex3u(1);
			}
		}
		else {
			int addp = 5;
			if (addflags[0]) {
				++addp;
				eBpi2(0,addp) = -ex3(0) * xyp2 * camerasnew[camerak].K(0,0) - ex3(1) * xyp2 * camerasnew[camerak].K(0,1);
				eBpi2(1,addp) = -ex3(1) * xyp2 * camerasnew[camerak].K(1,1);
			}
			if (addflags[1]) {
				++addp;
				eBpi2(0,addp) = -ex3(0) * xyp4 * camerasnew[camerak].K(0,0) - ex3(1) * xyp4 * camerasnew[camerak].K(0,1);
				eBpi2(1,addp) = -ex3(1) * xyp4 * camerasnew[camerak].K(1,1);
			}
			if (addflags[2]) {
				++addp;
				eBpi2(0,addp) = -ex3u(0);
				eBpi2(1,addp) = -ex3u(1);
			}
		}
		if (k == normcamera) {
			for(int j = 0; j < nrparir; ++j) {
				const double xp = eBpi2r(0,j);
				const double yp = eBpi2r(1,j);
				eBpi2rt(j,0) = xp * Wij[nro][0] + yp * Wij[nro][2];
				eBpi2rt(j,1) = xp * Wij[nro][2] + yp * Wij[nro][1];
			} // for
		} // if
		else {
			for(int j = 0; j < nrpari; ++j) {
				const double xp = eBpi2(0,j);
				const double yp = eBpi2(1,j);
				eBpi2t(j,0) = xp * Wij[nro][0] + yp * Wij[nro][2];
				eBpi2t(j,1) = xp * Wij[nro][2] + yp * Wij[nro][1];
			} // for
		} // else
		if (!solveflag) {
			if (k == normcamera) {
				eNpqqv[k][0].noalias() += eBpi2rt * eBpi2r;
				exiv[k].noalias() += eBpi2rt * efp2;
			}
			else {
				eNpqqv[k][0].noalias() += eBpi2t * eBpi2;
				exiv[k].noalias() += eBpi2t * efp2;
			}
		} // if !solveflag

		if (k == normcamera) {
			eNlij[ks] = eBpi2rt * eBpp;
		}
		else {
			eNlij[ks] = eBpi2t * eBpp;
		}

	} // if aoflag || k > 0
	if (!aoflag && k == 0 && nraddpi > 0) {
		int addp = -1;
		if (addflags[0]) {
			++addp;
			eBpi23(0,addp) = -ex3(0) * xyp2 * camerasnew[camerak].K(0,0) - ex3(1) * xyp2 * camerasnew[camerak].K(0,1);
			eBpi23(1,addp) = -ex3(1) * xyp2 * camerasnew[camerak].K(1,1);
		}
		if (addflags[1]) {
			++addp;
			eBpi23(0,addp) = -ex3(0) * xyp4 * camerasnew[camerak].K(0,0) - ex3(1) * xyp4 * camerasnew[camerak].K(0,1);
			eBpi23(1,addp) = -ex3(1) * xyp4 * camerasnew[camerak].K(1,1);
		}
		if (addflags[2]) {
			++addp;
			eBpi23(0,addp) = -ex3u(0);
			eBpi23(1,addp) = -ex3u(1);
		}
		for(int j = 0; j < nraddpi; ++j) {
			const double xp = eBpi23(0,j);
			const double yp = eBpi23(1,j);
			eBpit23(j,0) = xp * Wij[nro][0] + yp * Wij[nro][2];
			eBpit23(j,1) = xp * Wij[nro][2] + yp * Wij[nro][1];
		} // for

		if (!solveflag) {
			eNpqqv[0][0].noalias() += eBpit23 * eBpi23;
			exiv[0].noalias() += eBpit23 * efp2;
		} // if !solveflag

		eNlij[ks] = eBpit23 * eBpp;

	} // if k == 0

	for(int j = 0; j < 3; ++j) {
		const double xp = eBpp(0,j);
		const double yp = eBpp(1,j);
		// multiply the transposed part with the weight matrix
		eBppt(j,0) = xp * Wij[nro][0] + yp * Wij[nro][2];
		eBppt(j,1) = xp * Wij[nro][2] + yp * Wij[nro][1];
	} // for j

	eNppi.noalias() += eBppt * eBpp;
	etppi.noalias() += eBppt * efp2;

	return;
} // CalcDerivatives



void LSA::eCalcCalibratedDerivativesX(Eigen::MatrixXd &A, const int camerak) {

	const double xy2 = ex3(0) * ex3(0) + ex3(1) * ex3(1);
	const double xy4 = xy2 * xy2;
	const double k2 = camerasnew[camerak].distpara[0];
	const double k4 = camerasnew[camerak].distpara[1];
	const double xxpc = k4 * xy4 + ex3(0) * (2 * k2 * ex3(0) + 4 * k4 * ex3(0) * xy2) + k2 * xy2 + 1.;
	const double xypc = ex3(0) * (4 * k4 * xy2 * ex3(1) + 2 * k2 * ex3(1));
	const double yxpc = (4 * k4 * ex3(0) * xy2 + 2 * k2 * ex3(0)) * ex3(1);
	const double yypc = ex3(1) * (2 * k2 * ex3(1) + 4 * k4 * xy2 * ex3(1)) + k4 * xy4 + k2 * xy2 + 1.;

	int dim = (int) A.cols();
	if (dim == nrparir)
		dim = 5;
	if (dim == nrpari)
		dim = 6;
	for (int i = 0; i < dim; ++i) {
		const double xpo = A(0,i);
		const double ypo = A(1,i);

		const double xp = xxpc * xpo + xypc * ypo;
		const double yp = yypc * ypo + yxpc * xpo;
		A(0,i) = xp * camerasnew[camerak].K(0,0) + yp * camerasnew[camerak].K(0,1);
		A(1,i) = yp * camerasnew[camerak].K(1,1);
	}

	return;
} // eCalcCalibratedDerivativesX



void LSA::eCalcCalibratedDerivatives3(Eigen::Matrix<double, 2, 3> &A, const int camerak) {

	const double xy2 = ex3(0) * ex3(0) + ex3(1) * ex3(1);
	const double xy4 = xy2 * xy2;
	const double k2 = camerasnew[camerak].distpara[0];
	const double k4 = camerasnew[camerak].distpara[1];
	const double xxpc = k4 * xy4 + ex3(0) * (2 * k2 * ex3(0) + 4 * k4 * ex3(0) * xy2) + k2 * xy2 + 1.;
	const double xypc = ex3(0) * (4 * k4 * xy2 * ex3(1) + 2 * k2 * ex3(1));
	const double yxpc = (4 * k4 * ex3(0) * xy2 + 2 * k2 * ex3(0)) * ex3(1);
	const double yypc = ex3(1) * (2 * k2 * ex3(1) + 4 * k4 * xy2 * ex3(1)) + k4 * xy4 + k2 * xy2 + 1.;

	for (int i = 0; i < 3; ++i) {
		const double xpo = A(0,i);
		const double ypo = A(1,i);

		const double xp = xxpc * xpo + xypc * ypo;
		const double yp = yypc * ypo + yxpc * xpo;
		A(0,i) = xp * camerasnew[camerak].K(0,0) + yp * camerasnew[camerak].K(0,1);
		A(1,i) = yp * camerasnew[camerak].K(1,1);
	}

	return;
} // eCalcCalibratedDerivatives3



void LSA::CalcDerivatives3DGCP(const int i, const int gcpnum, const std::vector<double> &XM,
		const bool robustflagi, const bool itbreakflagi, const std::vector<std::vector<float> > &W3ijit,
		const std::vector<int> &gcptypeit, const std::vector<std::vector<double> > &gcpdatait,
		const bool derivflag) {

	double diff1, diff2, diff3;

	const int gpi = gcpinfo[i];

	if (itbreakflagi) {
		diff1 = gcpdatait[gpi][0] - XXnew(0,i);
		diff2 = gcpdatait[gpi][1] - XXnew(1,i);
		diff3 = gcpdatait[gpi][2] - XXnew(2,i);
	}
	else {
		diff1 = gcpdata[gpi][0] - XXnew(0,i);
		diff2 = gcpdata[gpi][1] - XXnew(1,i);
		diff3 = gcpdata[gpi][2] - XXnew(2,i);
	}

	sig0ao += (float) (diff1 * diff1 * W3ijit[0][gcpnum] + diff2 * diff2 * W3ijit[1][gcpnum] + diff3 * diff3 * W3ijit[2][gcpnum] +
			2 * (diff1 * diff2 * W3ijit[3][gcpnum] + diff1 * diff3 * W3ijit[4][gcpnum] + diff2 * diff3 * W3ijit[5][gcpnum]));
	if (robustflagi)
		sig0aorv[gcpnum] = (float) (diff1 * diff1 * W3ijold[0][gcpnum] + diff2 * diff2 * W3ijold[1][gcpnum] + diff3 * diff3 * W3ijold[2][gcpnum] +
				2 * (diff1 * diff2 * W3ijold[3][gcpnum] + diff1 * diff3 * W3ijold[4][gcpnum] + diff2 * diff3 * W3ijold[5][gcpnum]));

	if (!derivflag)
		return;

	Eigen::Vector3d efpp3(diff1, diff2, diff3);
	Eigen::Matrix3d eBpp3, eBpp3t;

	for (int j = 0; j < 3; ++j){
		for (int k = 0; k < 3; ++k)
			if (k == j)
				eBpp3(k,j) = -1.;
			else
				eBpp3(k,j) = 0.;
		const double xp = eBpp3(0,j);
		const double yp = eBpp3(1,j);
		const double zp = eBpp3(2,j);

		//multiply the transposed part with the weight matrix
		eBpp3t(j,0) = xp * W3ijit[0][gcpnum] + yp * W3ijit[3][gcpnum] + zp * W3ijit[4][gcpnum];
		eBpp3t(j,1) = xp * W3ijit[3][gcpnum] + yp * W3ijit[1][gcpnum] + zp * W3ijit[5][gcpnum];
		eBpp3t(j,2) = xp * W3ijit[4][gcpnum] + yp * W3ijit[5][gcpnum] + zp * W3ijit[2][gcpnum];
	} // for j

	eNppi.noalias() += eBpp3t * eBpp3;
	etppi.noalias() += eBpp3t * efpp3;

	return;
} // LSA::CalcDerivatives3DGCP



void LSA::DecodeCameras(const std::vector<double> &XM) {

	int nrimagesoffset = -1;
	if (aoflag)
		nrimagesoffset = 0;

	for (int k = 0; k < nrimages; ++k) {

		if (k != 0 || aoflag) {
			const int campk = camparv[k];
			const int kioff = k + nrimagesoffset;

			ComputeRotationMatrixfromRodrigues(XM[campk],XM[campk + 1],XM[campk + 2],eRabc);

			eRMnew[k] = eRabc * eRM[kioff];

			RectifyRotationMatrix3(eRMnew[k]);

			if (k == normcamera) { // normcamera
				eesd(esvm[0]) = XM[campk + 3];
				eesd(esvm[1]) = XM[campk + 4];
				eesd(esvmaxm) = 0.;
			} //if k == normcamera
			else
				for (int i = 0; i < 3; ++i)
					eesd(i) = XM[campk + 3 + i];
			ePCnew.col(k) = eesd + ePC.col(kioff);
			const int imgindk = imgind[k];
			if (k == normcamera) {
				int addp = 4;
				if (addflags[0])
					camerasnew[imgindk].distpara[0] = XM[++addp + campk];
				if (addflags[1])
					camerasnew[imgindk].distpara[1] = XM[++addp + campk];
				if (addflags[2])
					camerasnew[imgindk].K(0,0) = camerasnew[imgindk].K(1,1) = XM[++addp + campk];
			}
			else {
				int addp = 5;
				if (addflags[0])
					camerasnew[imgindk].distpara[0] = XM[++addp + campk];
				if (addflags[1])
					camerasnew[imgindk].distpara[1] = XM[++addp + campk];
				if (addflags[2])
					camerasnew[imgindk].K(0,0) = camerasnew[imgindk].K(1,1) = XM[++addp + campk];
			}
			camerasnew[imgindk].K(2,2) = 1.;
		} // if k
		else if (aoflag) {
			ePCnew.col(0) = ePCold0;
			eRMnew[0] = eRMold0;
		}
	} // for

	if (!aoflag) {
		const int imgind0 = imgind[0];
		int addp = -1;
		if (addflags[0])
			camerasnew[imgind0].distpara[0] = XM[++addp];
		if (addflags[1])
			camerasnew[imgind0].distpara[1] = XM[++addp];
		if (addflags[2])
			camerasnew[imgind0].K(0,0) = camerasnew[imgind0].K(1,1) = XM[++addp];
		camerasnew[imgind0].K(2,2) = 1.;
	}

	return;
} // LSA::DecodeCameras





// Function to solve equation system given in block form with Preconditioned Conjugate Gradients
bool LSA::SparsePCG(const int begini){

	const int nb = (int) eNpqqv.size();
	int msize = 0;

	if (iniflag) {
		eCv.resize(nb);
		erkv.resize(nb);
		erk1v.resize(nb);
		edkv.resize(nb);
		eAdkv.resize(nb);
		ehkv.resize(nb);
		ehk1v.resize(nb);
		ehkdv.resize(nb);
		exivold.resize(nb);
	}

	double erthk, edtAdk, er1th1, erth;

	for (int ii = 0; ii < nb; ++ii) {
		eCv[ii] = eNpqqv[ii][0];
		erkv[ii] = exiv[ii];
		erk1v[ii] = exiv[ii];
		edkv[ii] = exiv[ii];
		eAdkv[ii] = exiv[ii];
		ehkv[ii] = exiv[ii];
		ehk1v[ii] = exiv[ii];
		exivold[ii] = exiv[ii];
		exiv[ii].setZero(); // start with zero vector
	}

	// Generate block Jacobian for preconditioning and initialize conjugate gradients
	for (int i = begini; i < nb; ++i){

		eCv[i] = eNpqqv[i][0].inverse();

		msize += (int) exiv[i].rows();
		ehkv[i] = eCv[i] * erkv[i];
		edkv[i] = ehkv[i];
	}

	if (msize < 1000)
		msize = 1000;

	double meandrold = 1.e50, meandr0 = 1.e50;
	int k = 0;
	for (k = 0; k < msize; ++k) {

		eAAdv(begini, nb, eNpqqv, Npqqnumv, edkv, eAdkv);
		evvdv(begini, nb, erkv, ehkv, erthk);
		evvdv(begini, nb, edkv, eAdkv, edtAdk);
		double ealphak;
		if (fabs(erthk) < 1.e-20 || fabs(erthk) > 1.e20)
			ealphak = 0.;
		else
			ealphak = erthk / edtAdk;

		double meandr = 0.;
		for (int i = begini; i < nb; ++i){
			for (int j = 0; j < exiv[i].rows(); ++j) {
				exiv[i](j) += ealphak * edkv[i](j);
				const double dr = ealphak * eAdkv[i](j);
				erk1v[i](j) = erkv[i](j) - dr;
				meandr += dr * dr;
			}
			ehk1v[i] = eCv[i] * erk1v[i];
		}

		if (k == 0)
			meandr0 = meandr;

		if (meandr < meandrold) {
			if (exiv[1](0) > 1.e10 || exiv[1](0) < -1.e10)
				std::cout << "MM: " << k << ' ' << erthk << " " << ealphak << " " << meandr << " " << meandr0 << "   " << meandrold << " " << exiv.size() << '\n';
			for (int ii = 0; ii < nb; ++ii)
				exivold[ii] = exiv[ii];
			meandrold = meandr;
		}

		if (meandr * 1.e5 < meandr0)
			break;

		for (int i = begini; i < nb; ++i)
			ehkdv[i] = ehk1v[i] - ehkv[i];

		evvdv(begini, nb, erk1v, ehkdv, er1th1);
		evvdv(begini, nb, erkv, ehkv, erth);
		double ebetak;
		if (fabs(er1th1) < 1.e-10)
			ebetak = 0.;
		else
			ebetak = er1th1 / erth;
		if (ebetak < 0.)
			ebetak = 0.;

		for (int i = begini; i < nb; ++i){
			for (int j = 0; j < edkv[i].rows(); ++j)
				edkv[i](j) = ehk1v[i](j) + ebetak * edkv[i](j);
			erkv[i] = erk1v[i];
			ehkv[i] = ehk1v[i];
		}
	}

	for (int ii = 0; ii < nb; ++ii)
		exiv[ii] = exivold[ii];

	return (false);
} // SparsePCG



void LSA::PointXXvIn(const Eigen::MatrixXd &XXi, std::vector<std::vector<IMAPNR> > &pointXXvi,
		const std::vector<int> &gcpinfoin){
	int i;

	std::vector<IMPNR> pxv;

	pointXXv.resize(0);
	nrpv.resize(0);

	nrpts = 0;
	IMPNR impnr;

	gcpinfo.resize(0);

	std::vector<int> imgindbv((int) PMatrices.size(),-1);
	for (i = 0; i < (int) imgind.size(); ++i)
		imgindbv[imgind[i]] = i;

	std::vector<std::vector<IMAPNR> >::iterator pxxv;
	std::vector<IMAPNR>::iterator pxxiv;

	const bool gcpf = gcpflag || (int) gcpinfoin.size() > 0;

	for (pxxv = pointXXvi.begin(), i = 0; i < XXnr; ++pxxv, ++i) {
		int nriii = 0;
		pxv.resize(0);
		for (pxxiv = pxxv->begin(); pxxiv != pxxv->end(); ++pxxiv) {
			// iimage could be used (together with imgindv and nrpv etc. to
			// adjust only a given set of images
			if (pxxiv->image > -1 && imgindbv[pxxiv->image] > -1) {
				impnr.iimage = imgindbv[pxxiv->image];
				impnr.pointnr = pxxiv->pointnr;
				pxv.push_back(impnr);
				++nriii;
			}
		}
		if (fabs(XXi(2,i)) > 1.e-8 && (nriii >= minimagenr
				|| (gcpf && gcpinfoin[i] > 0 && nriii > 0))) {
			nrpv.push_back(i);
			pointXXv.push_back(pxv);
			if (gcpf)
				gcpinfo.push_back(gcpinfoin[i]);
			++nrpts;
		} // if nrii
	}

	XX = Eigen::MatrixXd(3,nrpts);
	XXnew = Eigen::MatrixXd(3,nrpts);

	for (i = 0; i < nrpts; ++i) {
		XX(0,i) = XXnew(0,i) = XXi(0,nrpv[i]);
		XX(1,i) = XXnew(1,i) = XXi(1,nrpv[i]);
		XX(2,i) = XXnew(2,i) = XXi(2,nrpv[i]);
	}

	return;
}



void LSA::PointXXvOut(Eigen::MatrixXd &XXi, std::vector<std::vector<IMAPNR> > &pointXXvio, std::vector<int> &gcpinfoo){
	std::vector<IMAPNR> pxvio;
	std::vector<IMPNR>::iterator pxxv;
	IMAPNR imapnr;

	const bool gcpf = gcpflag || (int) gcpinfo.size() > 0;

	int j = 0;
	for (int i = 0; i < XXnr; ++i) {
		if (j < (int) pointXXv.size() && nrpv[j] == i) {
			pxvio.resize(0);
			for (pxxv = pointXXv[j].begin(); pxxv != pointXXv[j].end(); ++pxxv) {
				imapnr.pointnr = pxxv->pointnr;
				imapnr.image = imgind[pxxv->iimage];
				pxvio.push_back(imapnr);
			} // for
			pointXXvio[i] = pxvio;
			XXi(0,i) = XX(0,j);
			XXi(1,i) = XX(1,j);
			XXi(2,i) = XX(2,j);
			if (gcpf)
				gcpinfoo[i] = gcpinfo[j];
			++j;
		} // if nrpv
		else {
			pointXXvio[i].clear();
			if (gcpf)
				gcpinfoo[i] = -1;
		}
	} // for i

	return;
}



void ComputeRotationMatrixfromRodrigues(const double a, const double b, const double c,
		Eigen::Matrix3d &Rabc){

	if (boost::math::isnan(a) || boost::math::isnan(b) || boost::math::isnan(c))
		std::cout << "NANRabc: " << a << ' ' << b << ' ' << c << '\n';

	const double a2 = a * a, b2 = b * b, c2 = c * c;
	const double ab = a * b, ac = a * c, bc = b * c;
	double abc = 4 + a2 + b2 + c2;
	if (fabs(abc) < 1.e-20) {
		if (abc < 0.)
			abc = -1.e20;
		else
			abc = 1.e20;
	}
	else
		abc = 1. / abc;

	Rabc(0,0) = abc * (4 + a2 - b2 - c2);
	Rabc(0,1) = abc * (2 * ab - 4 * c);
	Rabc(0,2) = abc * (2 * ac + 4 * b);
	Rabc(1,0) = abc * (2 * ab + 4 * c);
	Rabc(1,1) = abc * (4 - a2 + b2 - c2);
	Rabc(1,2) = abc * (2 * bc - 4 * a);
	Rabc(2,0) = abc * (2 * ac - 4 * b);
	Rabc(2,1) = abc * (2 * bc + 4 * a);
	Rabc(2,2) = abc * (4 - a2 - b2 + c2);

	return;
} // ComputeRotationMatrixfromRodrigues



// Block matrix multiplication for preconditioned conjugate gradients
inline void eAAdv(const int begini, const int nb, const std::vector<std::vector<Eigen::MatrixXd> > &A,
		const std::vector<std::vector<bool> > &Anum,
		const std::vector<Eigen::MatrixXd> &b,
		std::vector<Eigen::MatrixXd> &c) {

	for (int i = begini; i < nb; ++i)
		for (int j = 0; j < c[i].rows(); ++j)
			c[i](j) = 0.;

	for (int i = begini; i < nb; ++i) {
		for (int j = 0; j < (int) Anum[i].size(); ++j) {
			if (Anum[i][j]) {
				c[i-j].noalias() += A[i][j] * b[i];
				if (j > 0)
					c[i].noalias() += A[i][j].transpose() * b[i-j];
			} // if Anum
		} // for j
	} // for i

	return;
}



inline void evvdv(const int begini, const int nb, const std::vector<Eigen::MatrixXd> &a,
		const std::vector<Eigen::MatrixXd> &b,
		double &c){
	c = 0.;
	Eigen::Matrix<double, 1, 1> cc;

	for (int i = begini; i < nb; ++i) {
		cc = b[i].transpose() * a[i];
		c += cc(0);
	}

	return;
}
