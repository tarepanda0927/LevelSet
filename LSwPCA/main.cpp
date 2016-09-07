#include "naricommon.h"
#include "ip/narigaussian.h"
#include <iomanip>

#include "naritimer.h"
#include "naricaseinfo.h"
#include <algorithm>
#include "other/narimhd.h"

#include <sstream>
#include "naricommon.h"
#include "narisystem.h"
#include "naripathline.h"
#include "naricaseinfo.h"
#include "ip/narimorphology.h"
#include "ip/narilabeling.h"

#include "ip/naricontour.h"
#include "ip/naridistance.h"
#include "other/narimhd.h"
#include "naritimer.h"
#include "ip/narirbf.h"
#include "../mist/vector.h"
#include "narivectorpp.h"
#include "ip/nariinterpolate.h"
#include "info.h"
#include "raw_io.h"

int main(int argc, char *argv[])
{
	text_info input_info;
	input_info.input(argv[1]);
	//テキストデータ読み込み
	std::vector<std::string> fcase;
	std::vector<std::string> rcase;
	std::ifstream f_case(input_info.dir_list + input_info.case_flist);
	std::ifstream r_case(input_info.dir_list + input_info.case_rlist);
	std::string buf_ft;
	std::string buf_rt;
	while (f_case&& getline(f_case, buf_ft))
	{
		fcase.push_back(buf_ft);
	}
	while (r_case&& getline(r_case, buf_rt))
	{
		rcase.push_back(buf_rt);
	}

	//レベルセット
	//症例ループ
	for (int i = 0; i < fcase.size(); i++) {
		//ラベル画像読み込み
		nari::mhd mhdI;
		mhdI.load(input_info.dir_Fl + fcase[i] + "_label.mhd");
		int xe = mhdI.size1();
		int ye = mhdI.size2();
		int ze = mhdI.size3();
		nari::vector<unsigned char> imgF(xe * ye * ze), imgR(xe * ye * ze);
		imgF.load_file_bin(input_info.dir_Fl + fcase[i] + "_label.raw");
		imgR.load_file_bin(input_info.dir_Ref + rcase[i] + "_label.raw");

		//レベルセット画像を作成＆保存
		nari::vector<double> imgFD(xe * ye * ze), imgRD(xe * ye * ze);
		nari::distance::euclidean_signed_distance_transform(imgF.ptr(), imgFD.ptr(), xe, ye, ze);
		nari::distance::euclidean_signed_distance_transform(imgR.ptr(), imgRD.ptr(), xe, ye, ze);
		std::ostringstream ossf;
		std::ostringstream ossr;
		ossf << input_info.dir_out << "LS/Fl/";
		ossr << input_info.dir_out << "LS/Ref/";
		nari::system::make_directry(ossf.str());
		nari::system::make_directry(ossr.str());
		mhdI.save_mhd_and_image(imgFD, ossf.str() + fcase[i] + ".raw");
		mhdI.save_mhd_and_image(imgRD, ossr.str() + rcase[i] + ".raw");

		//3分割した画像を作成＆保存
		int L1 = input_info.a;
		int L2 = input_info.b;
		int z1 = L1 + 1;
		int z2 = L2 - L1;
		int z3 = ze - L2 + 1;

		nari::vector<double> imgFD1, imgFD2, imgFD3, imgRD1, imgRD2, imgRD3;
		for (int z = 0; z < ze; z++) {
			for (int y = 0; y < ye; y++) {
				for (int x = 0; x < xe; x++) {
					int s = xe*ye*z + xe*y + x;
					if (z < z1) {
						imgFD1.push_back(imgFD[s]);
						imgRD1.push_back(imgRD[s]);
					}
					if (z > z1 + 1 && z < z2){
						imgFD2.push_back(imgFD[s]);
						imgRD2.push_back(imgRD[s]);
					}
					if (z > z2 + 1) {
						imgFD3.push_back(imgFD[s]);
						imgRD3.push_back(imgRD[s]);
					}
				}
			}
		}
		imgFD1.save_file_bin(ossf.str() + fcase[i] + "_1.raw");
		imgRD1.save_file_bin(ossr.str() + rcase[i] + "_1.raw");
		imgFD2.save_file_bin(ossf.str() + fcase[i] + "_2.raw");
		imgRD2.save_file_bin(ossr.str() + rcase[i] + "_2.raw");
		imgFD3.save_file_bin(ossf.str() + fcase[i] + "_3.raw");
		imgRD3.save_file_bin(ossr.str() + rcase[i] + "_3.raw");
	}

	//wPCA
	//

}