/*�擪�œ��{���ł�����ł����΃\�[�X�c���[�ŕ\�������Ƃ��ɕ����������Ȃ��炵���̂�*/

#ifndef __INFO__
#define __INFO__

#include "nariinfocontroller.h"
#include "narifile.h"
#include "naricommon.h"
#include <string>

struct text_info
{
	std::string dir_Fl;
	std::string dir_Ref;
	std::string dir_out;
	std::string dir_list;
	std::string case_flist;
	std::string case_rlist;
	int a;
	int b;

	inline void input(const std::string &path) // �e�L�X�g������͏�����������
	{
		nari::infocontroller info;
		info.load(path);

		dir_Fl = nari::file::add_delim(info.get_as_str("dir_Fl"));
		dir_Ref = nari::file::add_delim(info.get_as_str("dir_Ref"));
		dir_out = nari::file::add_delim(info.get_as_str("dir_out"));
		dir_list = nari::file::add_delim(info.get_as_str("dir_txt"));
		case_flist = info.get_as_str("case_f");
		case_rlist = info.get_as_str("case_r");
		//3��������{�[�_�[���C����z�����W���(0�n�܂�)
		a = info.get_as_int("line_1");
		b = info.get_as_int("line_2");

		info.output();
	}
};


#endif