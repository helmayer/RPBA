/*

RPBA - Robust Parallel Bundle Adjustment

File rpbacore.cpp

Description: Parallel bundle adjustment core



Copyright 2019 Helmut Mayer, Bundeswehr University Munich, Germany, Helmut.Mayer@unibw.de

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/


#include "rpba.h"
#include "system.h"
#include <boost/math/special_functions/round.hpp>



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
		std::ostream& os, const int mt,
		const bool extevalflag) {

	std::vector<bool> addflags,addflagsold;
	addflags = addflagsin;

	int nraddp = 0;
	for (int i = 0; i < (int) addflags.size(); ++i)
		if (addflags[i])
			++nraddp;

	addflagsold = addflags;

	Eigen::Matrix3f Ksi;

	std::vector<std::vector<IMAPNR> > pointXXv,pointXXvold;
	pointXXv = pointXXvio;

	const int pxs = (int) pointXXv.size();
	int nrimages = (int) xg.size();

	bool gcpflag = false;

	Eigen::MatrixXd XXd(4,pxs), XXdold(4,pxs), XXdt(4,pxs), XXt(4,pxs);
	std::vector<Eigen::MatrixXd> Zs; // XXs;

	std::vector<int> imgind, imgindbv(pxs);
	std::vector<Camera> cameraps(nrimages);
	std::vector<int> gcpinfo;
	std::vector<std::vector<float> > W3iji(6);
	std::vector<PMatd> PMatricesd(nrimages),PMatricesdold(nrimages),PMatricesdt(nrimages);
	std::vector<std::vector<double> > gcpdata,gcpdatait;

	for (int i = 0; i < nrimages; ++i) {
		cameraps[i] = cameras[i];
		PMatricesd[i] = PMatrices[i].cast <double> ();
		PMatricesdold[i] = PMatrices[i].cast <double> ();
	}

	XXd = XX.cast <double> ();

	int minnrimages = 70;
	const int minnrthreads = boost::math::iround((float) nrimages / (float) minnrimages);

	int nrthreads = mt;

	if (minnrthreads < nrthreads)
		nrthreads = minnrthreads;

	float ihw2 = 0.f;
	for(int i = 0; i < nrimages; ++i)
		ihw2 += (float) ((iheight[i] + iwidth[i]) * (iheight[i] + iwidth[i])) * 0.25f;
	ihw2 /= (float) nrimages;

	bool parflag = nrthreads > 1 && nrimages > 5;

	if (extevalflag) {
		EvalLSA(nraddp, cameraps, XXd, xg, wx, wy, wxy, pointXXv, PMatricesd, ihw2);
		return;
	}


	if (parflag) {

		std::vector<float> xiout;
		Eigen::MatrixXf wap;
		std::vector<float> sig0vio;
		double sig0g = 0.;
		long redundancy = 0;

		gcpflag = true;

		std::cout << "#threads: " << nrthreads << " " << minnrthreads << '\n';

		std::vector<std::vector<bool> > pitflagsv(pxs);
		std::vector<LSA*> LSAv(nrthreads);
		std::vector<Eigen::MatrixXd> XXdd(nrthreads);
		std::vector<std::vector<PMatd> > PMatricesdd(nrthreads);
		std::vector<std::vector<Camera> > camerasoit(nrthreads);
		std::vector<std::vector<std::vector<IMAPNR> > > pointXXvv(nrthreads);
		std::vector<IMAPNR> pxvv;
		std::vector<std::vector<IMAPNR> > pxvvv(nrthreads);
		std::vector<int> nriv(nrthreads+1);
		std::vector<std::vector<int> > imgindv(nrthreads),gcpinfov(nrthreads),pindv(nrthreads);
		std::vector<std::vector<std::vector<std::vector<float> > > > w3dvivv(nrthreads);
		std::vector<std::vector<std::vector<float> > > sig0vvv(nrthreads);
		std::vector<std::vector<std::vector<float> > > W3ijiv(nrthreads);
		std::vector<std::vector<float> > sig0vvs(nrimages);
		std::vector<long> redundancyv(nrthreads);
		std::vector<float> sig0wv(nrimages);
		std::vector<double> sig0gv(nrthreads);
		std::vector<std::vector<std::vector<double> > > gcpdatav;
		std::vector<int> beginv(nrthreads),endv(nrthreads);
		std::vector<Camera> camerasnew(nrimages),camerasold(nrimages);


		std::vector<std::vector<int> > imagesperpart;
		std::vector<idx_t> part(nrimages);

		const System::Time startp = System::getTickCount();
		imagesperpart.resize(nrthreads);
		PartitionpointXXv(pointXXv, nrimages, nrthreads, part, imagesperpart);
		std::cout << "\n\nPartitioning runtime: " << (float) (System::getTickCount()-startp) / 1000.f << "s\n\n";

		nriv[0] = 0;
		for (int it = 0; it < nrthreads; ++it) {
			nriv[it+1] = nriv[it] + (int) imagesperpart[it].size();
			w3dvivv[it].resize(6);
		}
		nriv[nrthreads] = nrimages;

		std::vector<int> imageind(nrimages);

		int itj = -1;
		for (int it = 0; it < nrthreads; ++it)
			for (int j = 0; j < (int) imagesperpart[it].size(); ++j)
				imageind[++itj] = imagesperpart[it][j];

		for (int it = 0; it < nrthreads; ++it) {
			camerasoit[it] = cameras;
			PMatricesdd[it] = std::vector<PMatd>(nrimages);
			for (int i = 0; i < nrimages; ++i) {
				PMatricesdd[it][i] = PMatrices[i].cast <double> ();
			} // for i
			imgindv[it].resize(nriv[it+1] - nriv[it]);
			std::cout << "Thread: " << it << " #images " << nriv[it+1] - nriv[it] << " indices " << nriv[it] << ' ' << nriv[it+1] << '\n';
			for (int i = 0; i < nriv[it+1] - nriv[it]; ++i) {
				imgindv[it][i] = imageind[i + nriv[it]];
				imgindbv[imageind[i + nriv[it]]] = it;
			}
			std::sort(imgindv[it].begin(),imgindv[it].end());
			sig0vvv[it].resize(nrimages);
		} // for it


		gcpinfo.resize(pxs);
		std::vector<int>::iterator pxi;
		for (pxi = gcpinfo.begin(); pxi != gcpinfo.end(); ++pxi)
			(*pxi) = -1;

		sig0g = 0.;
		redundancy = 0;
		int nrparimages = (nrimages - 2) * 6 + 5 + nrimages * nraddp;


		int nrobs;
		double sig0opt = 1.e20;

		int iflag = 1, iter = 1;
		if (robustflag)
			iflag = 2;
		bool robustflagit;
		bool iniflag = true;

		gcpflag = true;

		for (int i = 0; i < pxs; ++i) {
			if ((int) pointXXv[i].size() > 0) {
				pitflagsv[i].resize(nrthreads);
				for (int j = 0; j < nrthreads; ++j)
					pitflagsv[i][j] = false;
				for (int j = 0; j < (int) pointXXv[i].size(); ++j)
					pitflagsv[i][imgindbv[pointXXv[i][j].image]] = true;
			} // if pointXXv
		} // for i

		gcpdatav.resize(nrthreads);

		float sig0frac = 1.01f;

		while (iflag > 0 && iter < 51) {
			robustflagit = false;

			if (iter > 1) {
				if (iter < 3)
					addflags[0] = addflags[1] = addflags[2] = false;
				else
					addflags = addflagsold;

				std::vector<std::ostringstream> osv(nrthreads);

				#pragma omp parallel for num_threads(nrthreads)
				for (int it = 0; it < nrthreads; ++it) {
					if (iter == 2)
						LSAv[it] = new LSA( camerasoit[it], addflags,
								true, imgindv[it], PMatricesdd[it],
								pointXXvv[it], XXdd[it], (int)pindv[it].size(), xg,
								wx, wy, wxy, W3iji,
								true, gcpinfov[it], gcpdatait,
								ihw2, false);
					else
						LSAv[it]->InputLSAP(addflags, camerasoit[it], imgindv[it], pointXXvv[it], PMatricesdd[it], XXdd[it], (int) pointXXvv[it].size(),
								gcpinfov[it]);

					LSAv[it]->Adjust(false, osv[it], true, 0, W3ijiv[it], true, gcpdatav[it]);

					LSAv[it]->OutputLSAP(PMatricesdd[it], pointXXvv[it], XXdd[it], gcpinfov[it], camerasoit[it]);
				} // for it

				for (int it = 0; it < nrthreads; ++it) {
					os << osv[it].str();
					for (int i = 0; i < (int) imgindv[it].size(); ++i) {
						PMatricesd[imgindv[it][i]] = PMatricesdd[it][imgindv[it][i]];
						const int iii = imgindv[it][i];
						cameraps[iii] = camerasoit[it][imgindv[it][i]];
						camerasnew[iii] = camerasoit[it][imgindv[it][i]];
					} // for i
				}
			} // it iter > 1

			if (iflag == 1) {
				robustflagit = robustflag;
				if (robustflag)
					sig0frac = 1.02f;
			}

			const int nrppt = pxs / nrthreads;
			if (iniflag) {
				for (int itt = 0; itt < nrthreads; ++itt) {
					beginv[itt] = nrppt * itt;
					endv[itt] = nrppt * (itt + 1);
					if (itt == nrthreads - 1)
						endv[itt] = pxs;

					const int eb = endv[itt] - beginv[itt];
					for (int vv = 0; vv < 6; ++vv) {
						w3dvivv[itt][vv].resize(eb);
						for (int ebb = 0; ebb < eb; ++ebb)
							w3dvivv[itt][vv][ebb].resize(nrthreads);
					}
					W3ijiv[itt].resize(6);

				} // for itt
			} // if iniflag

			if (iter > 1) {
				for (int i = 0; i < pxs; ++i)
					pointXXv[i].resize(0);

				for (int it= 0; it < nrthreads; ++it)
					for (int i = 0; i < (int) pindv[it].size(); ++i)
						for (size_t j = 0; j < pointXXvv[it][i].size(); ++j)
							pointXXv[pindv[it][i]].push_back(pointXXvv[it][i][j]);

				for (int i = 0; i < pxs; ++i)
					if ((int) pointXXv[i].size() < 2)
						pointXXv[i].clear();
			}
			else {
				for (int it = 0; it < nrthreads; ++it) {
					for (int i = 0; i < (int) imgindv[it].size(); ++i) {
						const int iii = imgindv[it][i];
						cameraps[iii] = cameras[iii];
						camerasnew[iii] =  cameras[iii];
					} // for i
				} // for it
			} // else

			for (int itt = 0; itt < nrthreads; ++itt)
				for (int i = 0; i < 6; ++i)
					W3ijiv[itt][i].resize(0);
			for (int i = 0; i < pxs; ++i) {
				if ((int) pointXXv[i].size() > 0) {
					for (int it = 0; it < nrthreads; ++it)
						pitflagsv[i][it] = false;
					for (int j = 0; j < (int) pointXXv[i].size(); ++j)
						pitflagsv[i][imgindbv[pointXXv[i][j].image]] = true;
				} // if pointXXv
			} // for i

			if (robustflagit) {
				#pragma omp parallel for num_threads(nrthreads)
				for (int itt = 0; itt < nrthreads; ++itt) {
					for (int iii = 0; iii < nrimages; ++iii)
						sig0vvv[itt][iii].resize(0);
					sig0gv[itt] = 0.;
					redundancyv[itt] = 0;
					for (int i = beginv[itt], ii = 0; i < endv[itt]; ++i, ++ii) {
						if ((int) pointXXv[i].size() > 0) {
							CalcLSAXMulti(cameraps,i,XXd,xg,wx,wy,wxy,pointXXv,PMatricesd,
									imgindbv, nrthreads, pitflagsv[i],
									w3dvivv[itt][0][ii],w3dvivv[itt][1][ii],w3dvivv[itt][2][ii],
									w3dvivv[itt][3][ii],w3dvivv[itt][4][ii],w3dvivv[itt][5][ii],
									sig0gv[itt],redundancyv[itt],sig0vvv[itt],sig0wv,1,false,false);
						} // if pointXX
					} // for i
				} // for itt

				sig0g = 0.;
				redundancy = 0;
				for (int itt = 0; itt < nrthreads; ++itt) {
					sig0g += sig0gv[itt];
					redundancy += redundancyv[itt];
				}

				for (int i = 0; i < nrimages; ++i)
					sig0vvs[i].resize(0);

				for (int itt = 0; itt < nrthreads; ++itt)
					for (int i = 0; i < nrimages; ++i)
						sig0vvs[i].insert(sig0vvs[i].begin(),sig0vvv[itt][i].begin(),sig0vvv[itt][i].end());
				for (int i = 0; i < nrimages; ++i) {
					if (sig0vvs[i].size() == 0) {
						std::cout << "Image " << i << " without observations!!!!!!!!!!!!!!!!!!!!!\n\n";
						sig0wv[i] = 0.f;
					}
					else {
						std::sort(sig0vvs[i].begin(),sig0vvs[i].end());
						sig0wv[i] = 1.f / (sig0vvs[i][(int) sig0vvs[i].size() / 2] * 2.1981f); // 1.4826 // 1.706f - 1.306^2
					}
				} // for i
			} // if robustflagit

			int robustmode = 0;
			if (robustflagit)
				robustmode = 2;


			bool weightonlyflag = false;
			if (iter == 1)
				weightonlyflag = true;

			#pragma omp parallel for num_threads(nrthreads)
			for (int itt = 0; itt < nrthreads; ++itt) {
				sig0gv[itt] = 0.;
				redundancyv[itt] = 0;
				for (int i = beginv[itt], ii = 0; i < endv[itt]; ++i, ++ii) {
					if ((int) pointXXv[i].size() > 0) {
						CalcLSAXMulti(cameraps,i,XXd,xg,wx,wy,wxy,pointXXv,PMatricesd,
								imgindbv, nrthreads, pitflagsv[i],
								w3dvivv[itt][0][ii],w3dvivv[itt][1][ii],w3dvivv[itt][2][ii],
								w3dvivv[itt][3][ii],w3dvivv[itt][4][ii],w3dvivv[itt][5][ii],
								sig0gv[itt],redundancyv[itt],sig0vvv[itt],sig0wv,robustmode,weightonlyflag,true);
					} // if pointxxv
				} // for i
			} // for itt


			sig0g = 0.;
			redundancy = 0;
			for (int itt = 0; itt < nrthreads; ++itt) {
				sig0g += sig0gv[itt];
				redundancy += redundancyv[itt];
			}

			nrobs = 0;
			for (int i = 0; i < pxs; ++i)
				nrobs += (int) pointXXv[i].size();
			nrobs *= 2;

			sig0g /= (double) (redundancy - nrparimages);

			if (iter > 0) {
				std::cout << "IT: " << iter << " " << iflag << "  " << sqrt(sig0g * ihw2) << ' ' << sig0opt / sig0g << "  ";
				std::cout << sig0g << " " << redundancy << " " << nrobs << " " << nrparimages << " " << nraddp << "  " << redundancy - nrparimages << '\n';
			}

			if (sig0g <= sig0opt) {
				XXdold = XXd;
				for (int i = 0; i < nrimages; ++i)
					PMatricesdold[i] = PMatricesd[i];
				pointXXvold = pointXXv;
				camerasold = camerasnew;
			}

			for (int itt = 0; itt < nrthreads; ++itt) {
				camerasoit[itt] = camerasnew;
				for (int i = 0; i < (int) imgindv[itt].size(); ++i)
					PMatricesdd[itt][imgindv[itt][i]] = PMatricesd[imgindv[itt][i]];
				pointXXvv[itt].resize(0);
				pindv[itt].resize(0);
			} // for itt

			int nrgcpnew = 0;
			for (int itt = 0; itt < nrthreads; ++itt)
				gcpdatav[itt].resize(0);

			for (int itt = 0; itt < nrthreads; ++itt) {
				for (int i = beginv[itt], ii = 0; i < endv[itt]; ++i, ++ii) {
					if ((int) pointXXv[i].size() < 2)
						continue;
					int gnum = 0;
					for (int it = 0; it < nrthreads; ++it) {
						if (pitflagsv[i][it])
							++gnum;
						if (gnum > 1)
							break;
					}
					if (gnum > 1) {
						gcpinfo[i] = nrgcpnew;
						for (int ittt = 0; ittt < nrthreads; ++ittt) {
							gcpdatav[ittt].push_back(std::vector<double>(3));
							gcpdatav[ittt][nrgcpnew][0] = XXd(0,i);
							gcpdatav[ittt][nrgcpnew][1] = XXd(1,i);
							gcpdatav[ittt][nrgcpnew][2] = XXd(2,i);
						} // for ittt
						for (int ittt = 0; ittt < nrthreads; ++ittt)
							if (pitflagsv[i][ittt]) {
								for (int vv = 0; vv < 6; ++vv)
									W3ijiv[ittt][vv].push_back(w3dvivv[itt][vv][ii][ittt]);
							}
						++nrgcpnew;
					}
					else {
						gcpinfo[i] = -1;
					} // else
				} // for i
			} // for itt


			for (int i = 0; i < pxs; ++i) {
				for (int it = 0; it < nrthreads; ++it) {
					pxvvv[it].resize(0);
				}
				for (int j = 0; j < (int) pointXXv[i].size(); ++j) {
					const int iipi = imgindbv[pointXXv[i][j].image];
					pxvvv[iipi].push_back(pointXXv[i][j]);
				}
				for (int it = 0; it < nrthreads; ++it) {
					if (pxvvv[it].size() > 0) {
						pointXXvv[it].push_back(pxvvv[it]);
						pindv[it].push_back(i);
					} // if
				} // for
			} // for i

			for (int it = 0; it < nrthreads; ++it)
				gcpinfov[it].resize((int) pindv[it].size());

			#pragma omp parallel for num_threads(nrthreads)
			for (int it = 0; it < nrthreads; ++it) {
				XXdd[it].resize(4,(int) pindv[it].size());
				for (int j = 0; j < (int) pindv[it].size(); ++j) {
					XXdd[it].col(j) = XXd.col(pindv[it][j]);
					gcpinfov[it][j] = gcpinfo[pindv[it][j]];
				} // for i
			} // for it

			if (iter > 1) {
				if (sig0g <= sig0opt) {
					if (sig0opt < sig0g * sig0frac && iter > 3) { // 1.05
						--iflag;
					}
					else {
						sig0opt = sig0g;
					}
				}
				else {
					if (iter > 3)
						--iflag;
				}
			}
			else
				iniflag = false;

			++iter;

		} // while iflag

		pointXXvio = pointXXvold;


		for (int i = 0; i < nrimages; ++i)
			PMatrices[i] = PMatricesdold[i].cast <float> ();

		XX = XXdold.cast <float> ();

		cameras = camerasold;

		for (int it = 0; it < nrthreads; ++it)
			delete LSAv[it];

	} // if nrimages >
	else {

		gcpflag = false;

		LSA *LSA1 = new LSA(cameras, addflags,
				false, imgind, PMatricesd, pointXXv, XXd, pxs, xg,
				wx, wy, wxy, W3iji,
				gcpflag, gcpinfo, gcpdata,
				ihw2, robustflag);

		LSA1->Adjust(robustflag, os, false, -1, W3iji, false, gcpdata);

		LSA1->OutputLSAc(PMatrices, pointXXvio, XX, cameras);

		delete LSA1;

	} // else


	return;
} //void CalcPLSAc



// Bundle adjustment for points for multiple images
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
		const std::vector<float> &sig0wv, const int robustmode, const bool weightonlyflag, const bool wxyflag) {

	int js;
	int iter = 1;
	int iflag = 1;
	if (robustmode > 0)
		iflag = 2;
	float sig0 = 1.e20f;
	float X0,X1,X2,Xold0,Xold1,Xold2,X0in,X1in,X2in;
	float sig0opt = 1.e20f;

	X0 = Xold0 = X0in = (float) XX(0,pi);
	X1 = Xold1 = X1in = (float) XX(1,pi);
	X2 = Xold2 = X2in = (float) XX(2,pi);

	bool convflag = false;

	const int pxs = (int) pointXXvio[pi].size();
	int pxsd = pxs;
	int pxs2 = pxs * 2;

	std::vector<float> p1v(pxs), p2v(pxs), p3v(pxs);
	std::vector<float> wxi(pxs),wyi(pxs),wxyi(pxs),wxio(pxs),wyio(pxs),wxyio(pxs),sig0v(pxs);

	Eigen::Matrix3d eNpp, eNppi;
	Eigen::Vector3d exip;
	Eigen::Vector2d ecorr2;
	Eigen::Matrix<double, 2, 3> eBpp2;
	Eigen::Matrix<double, 3, 2> eNpp12;
	Eigen::MatrixXd eBpp2wv(2 * pxs, 3);
	Eigen::MatrixXd eNppv;

	std::vector<IMAPNR>::iterator pxvi;
	std::vector<double> corrv(pxs2);
	std::vector<bool> delf(pxs,false);
	float wsum = 0.f;

	if (wxyflag)
		eNppv = Eigen::MatrixXd::Constant(nrthreads * 3, 3, -1.);

	for(pxvi = pointXXvio[pi].begin(), js = 0; pxvi != pointXXvio[pi].end(); ++pxvi, ++js){
		const int j = pxvi->image;
		const int pxx = pxvi->pointnr;
		wxi[js] = wxio[js] = wx[j][pxx];
		wyi[js] = wyio[js] = wy[j][pxx];
		wxyi[js] = wxyio[js] = wxy[j][pxx];

		wsum += wxi[js] + wyi[js];
	} // for pxvi

	CalcCorrPointMulti(cameraps, X0, X1, X2, pointXXvio[pi],corrv,PMatrices,xg,
			wxi,wyi,wxyi,wxio,wyio,wxyio,sig0opt,sig0v,p1v,p2v,p3v);

	if (weightonlyflag) {
		wsum /= (float) pxs2;
		sig0io += sig0opt / wsum;
		redundancy += pxs2 - 3;
	}

	while((iflag > -1) && (iter < 101)){ // 101

		eNpp.setZero();
		if (wxyflag)
			eNppv.setZero();
		for(pxvi = pointXXvio[pi].begin(), js = 0; pxvi != pointXXvio[pi].end(); ++pxvi, ++js) {
			const int j = pxvi->image;

			const double x = p1v[js];
			const double y = p2v[js];

			const double P3 = p3v[js];
			const double P1 = p1v[js] * p3v[js];
			const double P2 = p2v[js] * p3v[js];

			eBpp2(0,0) = PMatrices[j](0,0) * P3 - PMatrices[j](2,0) * P1;
			eBpp2(0,1) = PMatrices[j](0,1) * P3 - PMatrices[j](2,1) * P1;
			eBpp2(0,2) = PMatrices[j](0,2) * P3 - PMatrices[j](2,2) * P1;
			eBpp2(1,0) = PMatrices[j](1,0) * P3 - PMatrices[j](2,0) * P2;
			eBpp2(1,1) = PMatrices[j](1,1) * P3 - PMatrices[j](2,1) * P2;
			eBpp2(1,2) = PMatrices[j](1,2) * P3 - PMatrices[j](2,2) * P2;


			CalcCalibratedDerivatives(eBpp2(0,0), eBpp2(1,0), eBpp2(0,1), eBpp2(1,1), eBpp2(0,2), eBpp2(1,2),
					x, y, cameraps[j].distpara[0], cameraps[j].distpara[1],
					cameraps[j].K(0,0), cameraps[j].K(1,1), cameraps[j].K(0,1));

			const int js2 = 2 * js;
			const int js21 = js2 + 1;
			eBpp2wv(js2,0) = eBpp2(0,0) * wxi[js] + eBpp2(1,0) * wxyi[js];
			eBpp2wv(js21,0) = eBpp2(1,0) * wyi[js] + eBpp2(0,0) * wxyi[js];
			eBpp2wv(js2,1) = eBpp2(0,1) * wxi[js] + eBpp2(1,1) * wxyi[js];
			eBpp2wv(js21,1) = eBpp2(1,1) * wyi[js] + eBpp2(0,1) * wxyi[js];
			eBpp2wv(js2,2) = eBpp2(0,2) * wxi[js] + eBpp2(1,2) * wxyi[js];
			eBpp2wv(js21,2) = eBpp2(1,2) * wyi[js] + eBpp2(0,2) * wxyi[js];

			eNpp.noalias() += eBpp2.transpose() * eBpp2wv.block<2,3>(js2,0);

			if (wxyflag) {
				if (nrthreads == 1)
					continue;
				else { // if (nn0flag) {
					const int ij = imgindbv[j];
					for (int i = 0; i < nrthreads; ++i)
						if (pitflags[i] && i != ij) {
							eNppv.block<3,3>(3 * i,0).noalias() += eBpp2.transpose() * eBpp2wv.block<2,3>(js * 2,0);
						}
				} // else
			}

		}

		if (wxyflag) {
			if (nrthreads == 1) {
				eNppv = eNpp;
			}
		}

		const double m1 = 1.001;

		eNpp(0,0) *= m1;
		eNpp(1,1) *= m1;
		eNpp(2,2) *= m1;

		eNppi = eNpp.inverse();

		exip.setZero();
		for (int i = 0; i < pxs; ++i) {
			ecorr2(0) = corrv[i * 2];
			ecorr2(1) = corrv[i * 2 + 1];

			exip.noalias() += eNppi * eBpp2wv.block<2,3>(i * 2,0).transpose() * ecorr2;
		}

		X0 = Xold0 + (float) exip(0);
		X1 = Xold1 + (float) exip(1);
		X2 = Xold2 + (float) exip(2);

		CalcCorrPointMulti(cameraps, X0, X1, X2, pointXXvio[pi],corrv,PMatrices,xg,
				wxi,wyi,wxyi,wxio,wyio,wxyio,sig0,sig0v,p1v,p2v,p3v);

		if (iflag == 1 && robustmode == 2) {
			for(pxvi = pointXXvio[pi].begin(), js = 0; pxvi != pointXXvio[pi].end(); ++pxvi, ++js) {
				if (sig0v[js] * sig0wv[pxvi->image] > 16.f) {
					wxi[js] *= 1.e-4f;
					wyi[js] *= 1.e-4f;
					wxyi[js] *= 1.e-4f;
					if (wxi[js] + wyi[js] < 2.e-8f) {
						wxi[js] = wyi[js] = 1.e-8f;
						wxyi[js] = 0.f;
					}
					--pxsd;
					delf[js] = true;
				} // if sig0v
			} // for pxvi
			if (pxsd < 2)
				break;
		} // if robustflag

		if (sig0 < sig0opt * 1.001) {
			if (sig0opt < sig0 * 1.05f) // 1.01f
				--iflag;
			if (!weightonlyflag) {
				XX(0,pi) = Xold0 = X0;
				XX(1,pi) = Xold1 = X1;
				XX(2,pi) = Xold2 = X2;
			}
			sig0opt = sig0;
			convflag = true;
			if (wxyflag) {
				for (int i = 0; i < nrthreads; ++i) {
					const int i3 = 3 * i;
					const int i31 = i3 + 1;
					wXXv[i] = (float) eNppv(i3,0);
					wYYv[i] = (float) eNppv(i31,1);
					wZZv[i] = (float) eNppv(i3 + 2,2);
					wXYv[i] = (float) eNppv(i3,1);
					wXZv[i] = (float) eNppv(i3,2);
					wYZv[i] = (float) eNppv(i31,2);
				}
			}
		} //if sig0
		else {
			--iflag;
		}

		if (weightonlyflag)
			break;

		++iter;
	} //while iflag


	if (robustmode == 1) {
		for(pxvi = pointXXvio[pi].begin(), js = 0; pxvi != pointXXvio[pi].end(); ++pxvi, ++js)
			sig0vv[pxvi->image].push_back(sig0v[js]);
	} // if robustmode

	if (robustmode == 2) {
		std::vector<int> delv;
		for (int i = 0; i < pxs; ++i)
			if (delf[i])
				delv.push_back(i);
		if ((int) delv.size() > 0) {
			pxs2 -= (int) delv.size() * 2;
			for (int i = (int)(delv.size()) - 1; i > -1; --i) {
				pointXXvio[pi].erase(pointXXvio[pi].begin() + delv[i]);
			}
		}
	} // if robustmode

	if (pointXXvio[pi].size() < 2 || !convflag) {
		pointXXvio[pi].clear();
		convflag = false;
		XX(0,pi) = X0in;
		XX(1,pi) = X1in;
		XX(2,pi) = X2in;
	}
	else {
		if (!weightonlyflag) {
			wsum /= (float) pxs2;
			sig0io += sig0opt / wsum;
			redundancy += pxs2 - 3;
		}
	}

	return(convflag);
} // bool CalcLSAXMulti



// Computation of corrections for bundle adjustment for points only for multi image case
void CalcCorrPointMulti(const std::vector<Camera> &cameraps,
		const double X0, const double X1, const double X2,
		std::vector<IMAPNR> &pxv,
		std::vector<double> &corrv,
		const std::vector<PMatd> &PM,
		const std::vector<Eigen::MatrixXf> &xg,
		const std::vector<float> &wx, const std::vector<float> &wy, const std::vector<float> &wxy,
		const std::vector<float> &wxo, const std::vector<float> &wyo, const std::vector<float> &wxyo,
		float &sig0, std::vector<float> &sig0v,
		std::vector<float> &p1v, std::vector<float> &p2v, std::vector<float> &p3v){

	int js;

	sig0 = 0;
	std::vector<IMAPNR>::iterator pxvi;
	for(pxvi = pxv.begin(), js = 0; pxvi != pxv.end(); ++pxvi, ++js){
		const int j = pxvi->image;
		const int pxx = pxvi->pointnr;
		const int js2 = js * 2;

//		if (fabs(PM[j](0,3)) > 1.e5 || fabs(PM[j](1,3)) > 1.e5 || fabs(PM[j](2,3)) > 1.e5)
//			continue;

		double p1 = PM[j](0,0) * X0 + PM[j](0,1) * X1 + PM[j](0,2) * X2 + PM[j](0,3);
		double p2 = PM[j](1,0) * X0 + PM[j](1,1) * X1 + PM[j](1,2) * X2 + PM[j](1,3);
		double p3 = PM[j](2,0) * X0 + PM[j](2,1) * X1 + PM[j](2,2) * X2 + PM[j](2,3);

		if (fabs(p3) < 1.e-3) {
			if (p3 > 0.)
				p3 = 1.e3;
			else
				p3 = -1.e3;
		}
		else
			p3 = 1. / p3;

		p1 *= p3;
		p2 *= p3;

		p1v[js] = (float) p1;
		p2v[js] = (float) p2;
		p3v[js] = (float) p3;

		RadialReconstructionImage(p1, p2, cameraps[j]);

		p1 *= cameraps[j].K(0,0);
		p2 *= cameraps[j].K(0,0);

		const double dx = xg[j](0,pxx) - p1;
		const double dy = xg[j](1,pxx) - p2;

		corrv[js2] = dx;
		corrv[js2 + 1] = dy;
		const float sig0i = (float) (dx * dx * wx[js] + dy * dy * wy[js] + 2 * dx * dy * wxy[js]);
		sig0 += sig0i;
		sig0v[js] = (float) (dx * dx * wxo[js] + dy * dy * wyo[js] + 2 * dx * dy * wxyo[js]);
	}

	return;
} //void CalcCorrPointMulti



void CalcCalibratedDerivatives(double &xp0, double &yp0, double &xp1, double &yp1, double &xp2, double &yp2,
		const double x, const double y, const double k2, const double k4,
		const double fx, const double fy, const double s) {

	const double xpo0 = xp0;
	const double ypo0 = yp0;
	const double xpo1 = xp1;
	const double ypo1 = yp1;
	const double xpo2 = xp2;
	const double ypo2 = yp2;

	const double xy2 = x * x + y * y;
	const double xy4 = xy2 * xy2;

	const double xxpc = k4 * xy4 + x * (2 * k2 * x + 4 * k4 * x * xy2) + k2 * xy2 + 1.;
	const double xypc = x * (4 * k4 * xy2 * y + 2 * k2 * y);
	const double yxpc = y * (4 * k4 * x * xy2 + 2 * k2 * x);
	const double yypc = y * (2 * k2 * y + 4 * k4 * xy2 * y) + k4 * xy4 + k2 * xy2 + 1.;

	xp0 = (xxpc * xpo0 + xypc * ypo0);
	yp0 = (yypc * ypo0 + yxpc * xpo0);
	xp0 = xp0 * fx + yp0 * s;
	yp0 *= fy;

	xp1 = (xxpc * xpo1 + xypc * ypo1);
	yp1 = (yypc * ypo1 + yxpc * xpo1);
	xp1 = xp1 * fx + yp1 * s;
	yp1 *= fy;


	xp2 = (xxpc * xpo2 + xypc * ypo2);
	yp2 = (yypc * ypo2 + yxpc * xpo2);
	xp2 = xp2 * fx + yp2 * s;
	yp2 *= fy;

	return;
} // CalcCalibratedDerivatives



// Evaluation of bundle adjustment for calibrated images
void EvalLSA(const int nraddp, const std::vector<Camera> &cameraps,
		const Eigen::MatrixXd &XXd,
		const std::vector<Eigen::MatrixXf> &xg,
		const std::vector<std::vector<float> > &wxi, const std::vector<std::vector<float> > &wyi,
		const std::vector<std::vector<float> > &wxyi,
		std::vector<std::vector<IMAPNR> > &pointXXv,
		const std::vector<PMatd> &PMatricesd,
		const float ihw2) {

	long redundancy = 0;
	long nrobs = 0;
	float sig0g = 0.f;
	int ii,jj;
	const int nrimages =  (int) PMatricesd.size();
	std::vector<std::vector<IMAPNR> >::iterator pxv;
	std::vector<IMAPNR>::iterator pxvi;
	std::vector<float> p1v,p2v,p3v;

	for (pxv = pointXXv.begin(), ii = 0; pxv != pointXXv.end(); ++pxv, ++ii) {
		const int pxss = (int) pxv->size();
		const int pxss2 = pxss * 2;

		p1v.resize(pxss);
		p2v.resize(pxss);
		p3v.resize(pxss);

		if (pxss > 1) {
			std::vector<float> wxii(pxss),wyii(pxss),wxyii(pxss),sig0v(pxss);
			std::vector<double> corrv(pxss2);
			float wsum = 0.f;
			float sig0optf;
			int js,pxx;
			for(pxvi = pxv->begin(), js = 0; pxvi != pxv->end(); ++pxvi, ++js) {
				jj = pxvi->image;
				pxx = pxvi->pointnr;
				wxii[js] = wxi[jj][pxx];
				wyii[js] = wyi[jj][pxx];
				wxyii[js] = wxyi[jj][pxx];

				wsum += wxii[js] + wyii[js];
			} // for pxvi

			wsum = wsum / (float) pxss2;

			std::vector<bool> pitflags(1);
			std::vector<int> imgindbv;
			CalcCorrPointMulti(cameraps, XXd(0,ii), XXd(1,ii), XXd(2,ii), (*pxv),corrv,PMatricesd,xg,
					wxii,wyii,wxyii,wxii,wyii,wxyii,sig0optf,sig0v,p1v,p2v,p3v);
			sig0g += sig0optf / wsum;
			redundancy += pxss2 - 3;
			nrobs += pxss2;
		} // if pxss
	} // for pxv

	sig0g *= 1.f / (float) (redundancy - (nrimages - 2) * 6 - 5 - nrimages * nraddp);

	std::cout << "ITTTT: " << sqrtf(sig0g * ihw2) << ' ' << nrobs << " " << redundancy - (nrimages - 2) * 6 - 5 - nrimages * nraddp << '\n';

	return;
} // EvalLSA



void PartitionpointXXv(std::vector<std::vector<IMAPNR> > &pointXXv, const int nrimages,
		const int nrthreads, std::vector<idx_t> &part,
		std::vector<std::vector<int> > &imagesperpart)
{
  const int minpointnr = 0; // minimum number of points used for the graph


  const System::Time startp = System::getTickCount();


  std::vector<std::vector<int> > AM(nrimages);
  std::vector<int> BM(nrimages,0);
  for (int i = 0; i < nrimages; ++i) {
	  AM[i].resize(nrimages);
	  for (int j = 0; j < nrimages; ++j)
	  	  AM[i][j] = 0;
	  } // for i

  std::vector<std::vector<IMAPNR> >::iterator pxv;
  std::vector<IMAPNR>::iterator pxvi,pxvj;

  for (pxv = pointXXv.begin(); pxv != pointXXv.end(); ++pxv) {
	  if ((int) pxv->size() < 2)
		  continue;
	  for (pxvi = pxv->begin(); pxvi != pxv->end(); ++pxvi) {
		  const int fi = pxvi->image;
		  BM[fi] += (int) pxv->size();
		  for (pxvj = pxvi + 1; pxvj != pxv->end(); ++pxvj) {
			  const int si = pxvj->image;
			  ++AM[fi][si];
			  ++AM[si][fi];
		  }
	  }
  }

  int nredges = 0;
  for (int i = 0; i < nrimages; ++i) {
	  for (int j = 0; j < nrimages; ++j) {
		  if (AM[i][j] > minpointnr)
			  ++nredges;
	  }
  }

  std::vector<idx_t> xadj(nrimages + 1), adjncy(nredges), vwgt(nrimages);
  xadj[0] = 0;

  int indexadj = -1;
  for (int i = 0; i < nrimages; ++i) {

    int nrconn = 0;
    for (int j = 0; j < nrimages; ++j)
    	if (AM[i][j] > minpointnr) {
    		adjncy[++indexadj] = j;
    		++nrconn;
    	}

    xadj[i+1] = xadj[i] + nrconn;
    vwgt[i] = (idx_t) pow(BM[i],0.3333);
  } // for


  std::cout << "Preparation time: " << (float) (System::getTickCount()-startp) / 1000.f << "s\n";

  idx_t nrparts = nrthreads;
  idx_t inrimages = nrimages;
  idx_t edgecut;

  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_CONTIG] = 1;

  // One constraint
  idx_t ncon = 1;

  METIS_PartGraphRecursive(&inrimages, &ncon, &xadj[0], &adjncy[0],
		  &vwgt[0], nullptr, nullptr, &nrparts, nullptr,
		  nullptr, &options[0], &edgecut, &part[0]);


  std::cout << "After partitioning: " << (float) (System::getTickCount()-startp) / 1000.f << "s\n\n";

  for (int i = 0; i < nrimages; ++i)
	  imagesperpart[part[i]].push_back(i);

  std::cout << "Images per partition: ";
  for (int it = 0; it < nrthreads; ++it)
	  std::cout << imagesperpart[it].size() << ' ';
  std::cout << '\n';

  std::cout << "After analysis of parts: " << (float) (System::getTickCount()-startp) / 1000.f << "s\n";

  return;
} // PartitionpointXXv
