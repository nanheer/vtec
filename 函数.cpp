#include"����.h"
#include<math.h>
#include<iostream>
#include<fstream>
#include<vector>
//#include"matrix.h"
#include<iomanip>
#include<string>
#include <Eigen/Dense>
#include <io.h>
using Eigen::MatrixXd;
using namespace std;
/*#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
#define STD std
#else
#define STD
#endif

#ifndef _NO_TEMPLATE
typedef matrix<double> Matrix;
#else
typedef matrix Matrix;
#endif

#ifndef _NO_EXCEPTION
#  define TRYBEGIN()	try {
#  define CATCHERROR()	} catch (const STD::exception& e) { \
						cerr << "Error: " << e.what() << endl; }
#else
#  define TRYBEGIN()
#  define CATCHERROR()
#endif*/
//��ȡsp3�ļ�
void readsp3file(string strn, psp3 sp3file)
{
	s_sp3_ephe onesp3;
	ifstream nfile(strn, ios::in);
	if (!nfile){
		cout << "sp3�ļ��򿪴���" << endl;
		exit(0);
	}
	cout << "��ʼ��sp3�ļ�" << endl;
	for (int i = 0; i < 32; i++){//prn��ʼ��
		if (i < 9){
			sp3file->all_ephem[i].prn = "G0" + to_string(i + 1);
		}
		else{sp3file->all_ephem[i].prn = "G" + to_string(i + 1);}
	}
	string str1;
	getline(nfile, str1);
	while (str1.substr(0, 1) != "*")
	{
		getline(nfile, str1);
	}
	while (!nfile.eof())
	{
		if (str1.substr(0, 1) == "*"){
			onesp3.utime_n.year = atoi(str1.substr(3, 4).c_str());
			onesp3.utime_n.month = atoi(str1.substr(8, 2).c_str());
			onesp3.utime_n.day = atoi(str1.substr(11, 2).c_str());
			onesp3.utime_n.hour = atoi(str1.substr(14, 2).c_str());
			onesp3.utime_n.minute = atoi(str1.substr(17, 2).c_str());
			onesp3.utime_n.second = atof(str1.substr(20, 10).c_str());
		}
		else if(str1.substr(0, 2) == "PG"){
			onesp3.x = atof(str1.substr(5, 13).c_str())*1000.0;
			onesp3.y = atof(str1.substr(19, 13).c_str())*1000.0;
			onesp3.z = atof(str1.substr(33, 13).c_str())*1000.0;
			for (int i = 0; i < 32; i++){
				if (sp3file->all_ephem[i].prn == str1.substr(1, 3).c_str()){
					sp3file->all_ephem[i].sate_ephem.push_back(onesp3);
				}
			}
		}
		getline(nfile, str1);
	}
	nfile.close();
	cout << "sp3�ļ���ȡ���" << endl;
}

//��ȡBDGIM ģ�͵ķǷ���ϵ����Ԥ������
void read_a0(string path, double** ctable){
	ifstream ifile(path, ios::in);
	if (!ifile){
		cout << "a0�ļ��򿪴���" << endl;
		exit(0);
	}
	string str;
	int row = 0;
	getline(ifile, str);
	do{
		for (int i = 0; i < 17; i++){
			ctable[row][i] = atof(str.substr(1 + i * 16, 10).c_str())*pow(10,atof(str.substr(12 + i * 16, 4).c_str()));
		}
		row++;
		do{
			getline(ifile, str);
		} while (str.length() == 0 && !ifile.eof());
	} while (!ifile.eof());
	ifile.close();
	cout << "a0�ļ���ȡ���" << endl;
}
//��ȡo�ļ�
void readofile_vtec(string stro, pobs obsfile)//������λƽ��α����㷽ʽ��ȡ����ͬһ�����Ƿ�һ��
{
	cout << "���ļ���" << stro << endl;
	obs_value val;//һ��ʱ��һ�����ǵĹ۲�����
	od one_time_val;//һ��ʱ���������ǵĹ۲�����
	//GPS
	int pos[4] = { -1, -1, -1, -1 };//����������������е�λ��
	double value;//��Ӧ����ֵ
	int prn_num;//�еĲ�վ����������prn����G 1����ʽ,ͳһ��ΪG01����ʽ
	string  prn_name;
	ifstream ofile(stro, ios::in);
	if (!ofile){
		cout << stro<<"o�ļ��򿪴���" << endl;
		exit(0);
	}
	string str1, str_ins;//str_ins��ʱ���ݣ�����str2��str3��pushback����ʹ��
	vector<string> str2, str3;//str2�洢����prn��str3�洢һ�����ǵ���������
	getline(ofile, str1);
	//cout << str1 << endl;
	int a = 0;//���еĲ��������
	while (str1.substr(60, 13) != "END OF HEADER")
	{
		//cout << "��ʼ���ļ�ͷ"<< endl;
		//cout <<str1<< endl;
		//GPS
		if (str1.substr(60, 19) == "# / TYPES OF OBSERV")
		{
			a = atoi(str1.substr(4, 2).c_str());
			int b = 0;//��ֹһ������д�������еĲ��������
			//cout << "a��ֵ��" <<a<< endl;
			for (int i = 0; i < a; i++)
			{
				if (str1.substr(10 + 6 * (i - b), 2) == "C1") pos[0] = i;
				if (str1.substr(10 + 6 * (i - b), 2) == "P2") pos[1] = i;
				if (str1.substr(10 + 6 * (i - b), 2) == "L1") pos[2] = i;
				if (str1.substr(10 + 6 * (i - b), 2) == "L2") pos[3] = i;
				if ((i+1)%9 ==0) getline(ofile, str1), b = i+1;//һ�����д9�����͵�����
			}
			if (pos[0] == -1)cout << "û��C1����" << endl;
			if (pos[1] == -1)cout << "û��P2����" << endl;
			if (pos[2] == -1)cout << "û��L1����" << endl;
			if (pos[3] == -1)cout << "û��L2����" << endl;
		}

		if (str1.substr(60, 20) == "RINEX VERSION / TYPE"){
			obsfile->obsheaddata.obsdata_type = str1.substr(40, 1);
			//cout << str1.substr(60, 19) << endl;
		}
		if (str1.substr(60, 19) == "APPROX POSITION XYZ"){

			obsfile->obsheaddata.approx_coordinate.x = atof(str1.substr(1, 14).c_str());
			obsfile->obsheaddata.approx_coordinate.y = atof(str1.substr(15, 14).c_str());
			obsfile->obsheaddata.approx_coordinate.z = atof(str1.substr(29, 14).c_str());
			//cout << obsfile->obsheaddata.approx_coordinate.x << endl;
		}
		if (str1.substr(60, 20) == "ANTENNA: DELTA H/E/N"){
			obsfile->obsheaddata.antena_height.upping = atof(str1.substr(7, 7).c_str());
			obsfile->obsheaddata.antena_height.easting = atof(str1.substr(21, 7).c_str());
			obsfile->obsheaddata.antena_height.northing = atof(str1.substr(35, 7).c_str());
			//cout << obsfile->obsheaddata.antena_height.upping << endl;
		}
		if (str1.substr(60, 11) == "MARKER NAME"){
			obsfile->obsheaddata.station = str1.substr(0, 4);
			//cout << obsfile->obsheaddata.antena_height.upping << endl;
		}
		getline(ofile, str1);
	}
	//cout << "��ȡo�ļ�ͷ�ɹ�" << endl;
	//getline(ofile, str1);
	//int cou=0;
	while (!ofile.eof())
	{
		//getline(ofile, str1);
		// << "str1:" << str1 << endl;
		//cout << "str1.length():" << str1.length() << endl;
		//if (str1.empty() == true && ofile.eof()) break;//���һ��Ϊ������Ϊtrue���ַ���Ϊ��ʱsubstr�����ᱨ��
		do{
			getline(ofile, str1);
			//cout << "str1:" << str1 << endl;
		} while (str1.length() == 0 && !ofile.eof());//��ֹ��������һ�����߶������
		if (str1.empty() == true && ofile.eof()) break;//���һ��Ϊ������Ϊtrue���ַ���Ϊ��ʱsubstr�����ᱨ��
		//cout << "str1:" << str1 << endl;
		//if (str1.substr(1, 2) != "17")cout << "str1" << str1<<endl;
		if (str1.substr(1, 2) == "17" && (str1.substr(32, 1) == "G" || str1.substr(32, 1) == "R" || str1.substr(32, 1) == "E" || str1.substr(32, 1) == "C"))// == "17"��ζ�Ŷ�ȡ����17������ݣ����޸�
		//if (str1.substr(32, 1) == "G" || str1.substr(32, 1) == "R" || str1.substr(32, 1) == "E")
		{
			//cou++;
			//�����vector����
			str2.clear();
			one_time_val.prn_list.clear();
			one_time_val.one_obs_data.clear();
			int sat = 0;
			int b = 0;//��ֹһ������д�������еĲ��������
			one_time_val.utime_o.year = atoi(str1.substr(1, 2).c_str()) + 2000;
			one_time_val.utime_o.month = atoi(str1.substr(4, 2).c_str());
			one_time_val.utime_o.day = atoi(str1.substr(7, 2).c_str());
			one_time_val.utime_o.hour = atoi(str1.substr(10, 2).c_str());
			one_time_val.utime_o.minute = atoi(str1.substr(13, 2).c_str());
			one_time_val.utime_o.second = atof(str1.substr(16, 10).c_str());
			sat = atoi(str1.substr(30, 2).c_str());
			for (int k = 0; k < (sat-1) / 12; k++){//����12 ��������ʱ�����У��ٴ�һ�Ż���һ��sat/12�����ݣ�һ��12��
				getline(ofile, str_ins);
				str2.push_back(str_ins);
			}
			for (int i = 0; i < sat; i++)
			{
				//�����vector����
				str3.clear();
				val.p1 = 0.0; val.p2 = 0.0; val.l1 = 0.0; val.l2 = 0.0;
				if (i % 12 == 0 && i != 0){ b = i; str1 = str2[i / 12 - 1]; };//����
				//cout <<"i��ֵ" <<i<< endl;
				for (int k = 0; k < (a-1) / 5 + 1; k++){//ͬ��sat-1��һ�����ǵ��������͹۲����ݣ�һ��a/5�����ݣ�һ��5����������
					getline(ofile, str_ins);
					str3.push_back(str_ins);
				}
				if (str1.substr(32 + (i - b) * 3, 1) != "G") continue;//��֤��ȡ������ΪGPS
				prn_num = atoi(str1.substr(32 + (i - b) * 3+1, 2).c_str());//����С��10�����֣��еĲ�վ����������prn����G 1����ʽ,ͳһ��ΪG01����ʽ
				if (prn_num < 9){
					prn_name = "G0" + to_string(prn_num);
				}
				else{ prn_name = "G" + to_string(prn_num); }
				//cout << "prn:" << str1.substr(32 + (i - b) * 3, 3) << endl;
				for (int p = 0; p < 4; p++)
				{
					if (str3[pos[p] / 5].length() < (14 + 16 * (pos[p] % 5)))break;//��֤��ȡ�����ݳ��ȹ�������Ϊatof�����Ĵ���ֵ����Ϊ��
					value = atof(str3[pos[p] / 5].substr(1 + 16 * (pos[p] % 5), 14).c_str());//pos[p]������������������Ϊpos[p]/5������Ϊpos[p]%5
					switch (p)
					{
					case 0:val.p1 = value; break;
					case 1:val.p2 = value; break;
					case 2:val.l1 = value; break;
					case 3:val.l2 = value; break;
					}
				}
				if (val.p1 == 0 || val.p1 < 20000000.0 || val.p2 == 0 || val.p2 < 20000000.0 || val.l1 == 0 || val.l2 == 0) continue;//�����˿�����
				one_time_val.prn_list.push_back(prn_name);
				one_time_val.one_obs_data.push_back(val);
			}
			obsfile->obsdata.push_back(one_time_val);
		}
		//getline(ofile, str1);//�������һ�з�ֹ���ֶ�ȡ������substr�����жϳ���
	}
	ofile.close();
	cout << "��ȡo�ļ��ɹ� " << endl;
	if (obsfile->obsdata.size() != 2880){ cout << "                                             ��Ԫ��������2880   " << obsfile->obsdata.size() << endl; }
}
/*void readofile_vtec(string stro, pobs obsfile)//������λƽ��α����㷽ʽ��ȡ����ͬһ�����Ƿ�һ��
{
	for (int i = 0; i < 32; i++){//prn��ʼ��
		if (i < 9){
			obsfile->all_sate[i].prn = "G0" + to_string(i + 1);
		}
		else{ obsfile->all_sate[i].prn = "G" + to_string(i + 1); }
	}
	obs_value val;
	//GPS
	int pos[4] = { -1,-1,-1,-1};//����������������е�λ��
	double value;//��Ӧ����ֵ
	
	ifstream ofile(stro, ios::in);
	if (!ofile){
		cout << "o�ļ��򿪴���" << endl;
		exit(0);
	}
	//memset(&obsfile->obsdata, 0, sizeof(obsfile->obsdata));
	//cout << "��ʼ��o�ļ�" << endl;
	string str1,str_ins;//str_ins��ʱ���ݣ�����str2��str3��pushback����ʹ��
	vector<string> str2, str3;//str2�洢����prn��str3�洢һ�����ǵ���������
	getline(ofile, str1);
	//cout << str1 << endl;
	int a = 0;//���еĲ��������
	while (str1.substr(60, 13) != "END OF HEADER")
	{
		//cout << "��ʼ���ļ�ͷ"<< endl;
		//cout <<str1<< endl;
		//GPS
		if (str1.substr(60, 19) == "# / TYPES OF OBSERV")
		{
			a = atoi(str1.substr(4, 2).c_str());
			int b = 0;//��ֹһ������д�������еĲ��������
			//cout << "a��ֵ��" <<a<< endl;
			for (int i = 0; i < a; i++)
			{
				if (str1.substr(10 + 6 * (i - b), 2) == "C1") pos[0] = i;
				if (str1.substr(10 + 6 * (i - b), 2) == "P2") pos[1] = i;
				if (str1.substr(10 + 6 * (i - b), 2) == "L1") pos[2] = i;
				if (str1.substr(10 + 6 * (i - b), 2) == "L2") pos[3] = i;
				if (i == 8) getline(ofile, str1), b = 9;//һ�����д9�����͵�����
			}

			if (pos[0] == -1)cout << "û��C1����" << endl;
			if (pos[1] == -1)cout << "û��P2����" << endl;
			if (pos[2] == -1)cout << "û��L1����" << endl;
			if (pos[3] == -1)cout << "û��L2����" << endl;
		}

		if (str1.substr(60, 20) == "RINEX VERSION / TYPE"){
			obsfile->obsheaddata.obsdata_type = str1.substr(40, 1);
			//cout << str1.substr(60, 19) << endl;
		}
		if (str1.substr(60, 19) == "APPROX POSITION XYZ"){

			obsfile->obsheaddata.approx_coordinate.x = atof(str1.substr(1, 14).c_str());
			obsfile->obsheaddata.approx_coordinate.y = atof(str1.substr(15, 14).c_str());
			obsfile->obsheaddata.approx_coordinate.z = atof(str1.substr(29, 14).c_str());
			//cout << obsfile->obsheaddata.approx_coordinate.x << endl;
		}
		if (str1.substr(60, 20) == "ANTENNA: DELTA H/E/N"){
			obsfile->obsheaddata.antena_height.upping = atof(str1.substr(7, 7).c_str());
			obsfile->obsheaddata.antena_height.easting = atof(str1.substr(21, 7).c_str());
			obsfile->obsheaddata.antena_height.northing = atof(str1.substr(35, 7).c_str());
			//cout << obsfile->obsheaddata.antena_height.upping << endl;
		}
		if (str1.substr(60, 11) == "MARKER NAME"){
			obsfile->obsheaddata.station = str1.substr(0, 4);
			//cout << obsfile->obsheaddata.antena_height.upping << endl;
		}
		getline(ofile, str1);
	}
	cout << "��ȡo�ļ�ͷ�ɹ�" << endl;
	//getline(ofile, str1);
	while (!ofile.eof())
	{
		//getline(ofile, str1);
		// << "str1:" << str1 << endl;
		//cout << "str1.length():" << str1.length() << endl;
		//if (str1.empty() == true && ofile.eof()) break;//���һ��Ϊ������Ϊtrue���ַ���Ϊ��ʱsubstr�����ᱨ��
		do{
			getline(ofile, str1);
			//cout << "str1:" << str1 << endl;
		} while (str1.length() == 0 && !ofile.eof());//��ֹ��������һ�����߶������
		if (str1.empty() == true && ofile.eof()) break;//���һ��Ϊ������Ϊtrue���ַ���Ϊ��ʱsubstr�����ᱨ��
		//cout << "str1:" << str1 << endl;
		if (str1.substr(1, 2) == "17"&&str1.length()>=58)// == "17"��ζ�Ŷ�ȡ����17������ݣ����޸�
		{
			//�����vector����
			str2.clear();
			int sat = 0;
			int b = 0;//��ֹһ������д�������еĲ��������
			val.utime_o.year = atoi(str1.substr(1, 2).c_str())+2000;
			val.utime_o.month = atoi(str1.substr(4, 2).c_str());
			val.utime_o.day = atoi(str1.substr(7, 2).c_str());
			val.utime_o.hour = atoi(str1.substr(10, 2).c_str());
			val.utime_o.minute = atoi(str1.substr(13, 2).c_str());
			val.utime_o.second = atof(str1.substr(16, 10).c_str());
			sat = atoi(str1.substr(30, 2).c_str());
			for (int k = 0; k < sat / 12; k++){//һ��sat/12�����ݣ�һ��12��
				getline(ofile, str_ins);
				str2.push_back(str_ins);
			}
			//cout << sat / 12 + 1 << endl;
			//if (sat > 12) getline(ofile, str2);//һ������д�������еĹ۲⵽������
			//cout << obsd.obs_n << endl;
			//cout << "�۲�ʱ�䣺" << val.utime_o.year << ":" << val.utime_o.month << ":" << val.utime_o.day << ":" << val.utime_o.hour << ":" << val.utime_o.minute << ":" << val.utime_o.second << endl;
			for (int i = 0; i < sat; i++)
			{
				//�����vector����
				str3.clear();
				val.p1 = 0.0; val.p2 = 0.0; val.l1 = 0.0; val.l2 = 0.0;
				if (i % 12 == 0 && i != 0){ b = i; str1 = str2[i / 12 - 1];};
				//cout <<"i��ֵ" <<i<< endl;
				for (int k = 0; k < a / 5+1; k++){//һ��a/5�����ݣ�һ��5����������
					getline(ofile, str_ins);
					str3.push_back(str_ins);
				}
				cout << "str3:" << str3 << endl;
				cout << "str4:" << str4 << endl;
				cout << "str3.length():" << str3.length() << endl;
				cout << "str4.length():" << str4.length() << endl;
				//GPS
				cout << "(i+1)%12:" << (i + 1) % 12 << endl;
				cout << "b:" << b<< endl;
				cout << "32 + (i - b) * 3:" << 32 + (i - b) * 3 << endl;
				cout << "str1:" << str1 << endl;
				cout << "str1.substr(32 + (i - b) * 3, 1):" << str1.substr(32 + (i - b) * 3, 1) << endl;
				//cout << "str1:" << str1 << endl;
				if (str1.substr(32 + (i - b) * 3, 1) != "G" ) continue;//��֤��ȡ������ΪGPS
				//cout << "prn:" << str1.substr(32 + (i - b) * 3, 3) << endl;
				for (int p = 0; p < 4; p++)
				{
					if (str3[pos[p]/5].length() < (14 + 16 * (pos[p]%5)))break;//��֤��ȡ�����ݳ��ȹ�������Ϊatof�����Ĵ���ֵ����Ϊ��
					value = atof(str3[pos[p] / 5].substr(1 + 16 * (pos[p]%5), 14).c_str());//pos[p]������������������Ϊpos[p]/5������Ϊpos[p]%5
					switch (p)
					{
					case 0:val.p1 = value; break;
					case 1:val.p2 = value; break;
					case 2:val.l1 = value; break;
					case 3:val.l2 = value; break;
					}
				}
				//cout << "hhh" << endl;
				if (val.utime_o.hour == 9 && val.utime_o.minute == 26 && str1.substr(32 + (i - b) * 3, 3) == "G07"){
					cout.setf(ios::fixed);
					cout << "val.p1:" << val.p1 << endl;
					cout << "val.p2:" << val.p2 << endl;
					cout << "val.l1:" << val.l1 << endl;
					cout << "val.l2:" << val.l2 << endl;
				}
				if (val.p1 == 0 || val.p1 < 20000000.0 || val.p2 == 0 || val.p2 < 20000000.0 || val.l1 == 0 || val.l2 == 0) continue;
				//if (val.p1 == 0 || val.p1 < 20000000.0 || val.p2 == 0 || val.p2 < 20000000.0 || val.l1 == 0 || val.l1*c / f1 < 20000000.0 || val.l2 == 0 || val.l2*c / f2 < 20000000.0) continue;//==0ʱ�ַ���Ϊ�գ���ֹΪ���ַ�����< 20000000.0֤����λֵ�����⣬���ǹ�������߶��Ǵ��������׵�
				//if (val.p1 < 20000000.0 || val.p2 < 20000000.0 ||  val.l1*c / f1 < 20000000.0  || val.l2*c / f2 < 20000000.0)cout << "����̫С" << endl;
				for (int j = 0; j < 32; j++)
				{
					//cout << "prn:" << str1.substr(32 + (i - b) * 3, 3) << endl;
					//cout << "obsfile->all_sate[j].prn:" << obsfile->all_sate[j].prn << endl;
					if (atoi(obsfile->all_sate[j].prn.substr(1, 2).c_str()) == atoi(str1.substr(32 + (i - b) * 3+1, 2).c_str()))//�еĲ�վ����������prn����G 1����ʽ
					{
						(obsfile->all_sate[j].obs_data.push_back(val));
					}
				}
			}
		}
		//getline(ofile, str1);//�������һ�з�ֹ���ֶ�ȡ������substr�����жϳ���
	}
	ofile.close();
	cout << "��ȡo�ļ��ɹ� " << endl;
}*/
//��ȡ���������ļ�
void read_ionex(string stri, pio ionfile)//������λƽ��α����㷽ʽ��ȡ����ͬһ�����Ƿ�һ��
{
	cout << "��ʼ��ionex�ļ� " << endl;
	satedcb dcb_sate;
	stadcb dcb_sta;
	vtecmap vtec_val;
	//vtecmap val_rms;
	ifstream ofile(stri, ios::in);
	int lat = 0;//γ�ȸ���
	int num;//һ�����ݰ����ĵ�����������9��������16������73��
	if (!ofile){
		cout << "ionex�ļ��򿪴���" << endl;
		exit(0);
	}
	//memset(&obsfile->obsdata, 0, sizeof(obsfile->obsdata));
	//cout << "��ʼ��o�ļ�" << endl;
	string str1, str2;
	getline(ofile, str1);
	//cout << str1 << endl;
	while (str1.substr(60, 13) != "END OF HEADER")
	{
		if (str1.substr(60, 16) == "PRN / BIAS / RMS")
		{
			//cout << "����dcb��" << str1 << endl;
			dcb_sate.prn = str1.substr(3, 3);
			dcb_sate.bias = atof(str1.substr(10, 6).c_str());
			dcb_sate.rms = atof(str1.substr(20, 6).c_str());
			ionfile->sate_dcb.push_back(dcb_sate);
		}
		if (str1.substr(60, 20) == "STATION / BIAS / RMS")
		{
			dcb_sta.station = str1.substr(6, 4);
			dcb_sta.bias = atof(str1.substr(30, 6).c_str());
			dcb_sta.rms = atof(str1.substr(40, 6).c_str());
			ionfile->sta_dcb.push_back(dcb_sta);
		}
		getline(ofile, str1);
	}
	cout << "��ȡionex�ļ�ͷ�ɹ�" << endl;
	//cout <<str1<< endl;
	while (str1.substr(60, 11) != "END OF FILE")//getline�����������в������¶�
	{
		getline(ofile, str1);//����END OF TEC MAP����END OF RMS MAP
		lat = 0;
		if (str1.substr(60, 16) == "START OF TEC MAP" || str1.substr(60, 16) == "START OF RMS MAP")//�µ�ʱ����Ӧ��ȫ��vtec
		{
			str2 = str1;//�����ж���START OF TEC MAP����START OF RMS MAP
			getline(ofile, str1);
			if (str1.substr(60, 20) == "EPOCH OF CURRENT MAP"){
				//cout << str1 << endl;
				vtec_val.vtime.year = atoi(str1.substr(2, 4).c_str());
				vtec_val.vtime.month = atoi(str1.substr(10, 2).c_str());
				vtec_val.vtime.day = atoi(str1.substr(16, 2).c_str());
				vtec_val.vtime.hour = atoi(str1.substr(22, 2).c_str());
				vtec_val.vtime.minute = atoi(str1.substr(28, 2).c_str());
				vtec_val.vtime.second = atof(str1.substr(34, 10).c_str());
			}
			getline(ofile, str1);
			while (str1.substr(60, 20) == "LAT/LON1/LON2/DLON/H")//�µ�γ�ȶ�Ӧ��vtec
			{
				num = 16;
				for (int i = 0; i < 5; i++)
				{
					if (i == 4) num = 9;
					getline(ofile, str1);
					for (int j = 0; j < num; j++)
					{
						vtec_val.tec_values[lat][i * 16 + j] = atof(str1.substr(j * 5, 5).c_str());
					}
				}
				lat++;
				getline(ofile, str1);
			}
			if (str2.substr(60, 16) == "START OF TEC MAP"){
				ionfile->vtec_map.push_back(vtec_val);//vtecֵ
			}
			else{
				ionfile->vtec_rms.push_back(vtec_val);//vtec���ֵ
			}
			//getline(ofile, str1);//����END OF TEC MAP����END OF RMS MAP
		}
	}
	ofile.close();
	cout << "��ȡionex�ļ��ɹ� " << endl;
}
//��ȡ�������г�����ļ�
//��ȡ�ļ�
void read_sh(string path, psh sh_file)
{
	ifstream  ifile(path, ios::in);
	one_coeff coe;
	one_line_coeff one_coe;
	if (!ifile){
		cout << "������ļ��򿪴���" << endl;
		exit(0);
	}
	string str;
	getline(ifile, str);
	do{
		//����ͷ
		while (str.substr(0, 42) != "DEGREE  ORDER    VALUE (TECU)   RMS (TECU)")
		{
			if ((str.substr(0, 41) == "COORDINATES OF EARTH-CENTERED DIPOLE AXIS"))
			{
				getline(ifile, str);
				coe.pole_lat = atof(str.substr(50, 6).c_str());
				getline(ifile, str);
				coe.pole_lon= atof(str.substr(50, 6).c_str());
			}
			if ((str.substr(0, 18) == "PERIOD OF VALIDITY"))
			{
				getline(ifile, str);
				coe.jul.year = atoi(str.substr(49, 4).c_str());
				coe.jul.month = atoi(str.substr(54, 2).c_str());
				coe.jul.day = atoi(str.substr(57, 2).c_str());
				coe.jul.hour = atoi(str.substr(60, 2).c_str());
				coe.jul.minute = atoi(str.substr(63, 2).c_str());
				coe.jul.second = atof(str.substr(66, 2).c_str());
			}
			getline(ifile, str);
			if (str.empty() == true && ifile.eof()) break;
		}
		getline(ifile, str);
		coe.one_coe.clear();
		//������
		while (str.empty() != true)//ÿ�θ�������������һ������
		{
			one_coe.degree = atoi(str.substr(2, 2).c_str());
			one_coe.order = atoi(str.substr(9, 3).c_str());
			one_coe.tec = atof(str.substr(18, 12).c_str());
			coe.one_coe.push_back(one_coe);
			getline(ifile, str);
		}
		sh_file->coeff.push_back(coe);
		//ÿ���������һ��
		while (str.empty() == true){//��ֹ�����ֶ������
			getline(ifile, str);
			if (ifile.eof()) break;
		}
	} while (!ifile.eof());
	ifile.close();
	cout << "��г����ϵ����ȡ���" <<endl;
}
//��ȡdcb�ļ�
void read_dcb(string strd, pdcb dcbfile){
	cout << "��ʼ��dcb�ļ� " << endl;
	dcb_f one_dcb;//һ������
	ifstream ofile(strd, ios::in);
	if (!ofile){
		cout << "dcb�ļ��򿪴���" << endl;
		exit(0);
	}
	string str1;
	int n = 0;
	//cout << str1 << endl;
	while (!ofile.eof())
	{
		getline(ofile, str1);
		//cout << "dcb�ļ�������" << str1.length() << endl;
		//cout << "str1.substr(1, 3)��" << str1.substr(1, 3) << endl;
		//cout << str1 << endl;
		if ((str1.substr(1, 3) == "DSB" || str1.substr(1, 3) == "DCB") && str1.length()>100)
		{
			//cout << "dcb�ļ�������" << ++n << endl;
			one_dcb.prn = str1.substr(11, 3);
			one_dcb.station = str1.substr(15, 4);
			//cout << "str1.substr(11, 3):" << one_dcb.prn << "  str1.substr(15, 4):" << one_dcb.station << endl;
			one_dcb.obs1 = str1.substr(25, 3);
			one_dcb.obs2 = str1.substr(30, 3);
			one_dcb.bias = atof(str1.substr(82, 9).c_str());
			one_dcb.rms = atof(str1.substr(94, 9).c_str());
			dcbfile->dcb_val.push_back(one_dcb);
		}
		if (str1.substr(0, 8) == "%=ENDBIA")break;
	}
	ofile.close();
	cout << "��ȡdcb�ļ��ɹ� " << endl;
	cout << "dcb�ļ����������" << dcbfile->dcb_val.size()<< endl;
	//for (int i = 0; i < dcbfile->dcb_val.size(); i++){
	//	cout << "prn" << dcbfile->dcb_val[i].prn << "station" << dcbfile->dcb_val[i].station << "obs1" << dcbfile->dcb_val[i].obs1 << "obs2" << dcbfile->dcb_val[i].obs2 << "dcb" << dcbfile->dcb_val[i].bias << endl;
//	}
}
void gpsttoutc(pgpst gt, ptc ut)
{
	pjulian ju;
	ju = (pjulian)malloc(sizeof(JULIANDAY));
	gpsttojulianday(gt,ju);
	juliandaytoutc(ju, ut);
	free(ju);
}
void juliandaytobdt(pjulian ju, pbdt bt)
{
	double dt = ju->daynum + (ju->secondfrac + ju->secondnum) / _DAY_IN_SECOND - 2453736.500000;
	bt->week = int(dt / 7);
	bt->second = (dt - bt->week*7)*_DAY_IN_SECOND;
}
void utctobdt(ptc ut, pbdt bt)
{
	if (ut->year < 2006 || ut->month>12 || ut->month < 0 || ut->day>31 || ut->day < 0 || ut->hour>24 || ut->hour < 0 || ut->minute>60 || ut->minute < 0 || ut->second>60 || ut->second < 0)
	{
		cout << "ʱ���С����"<<endl;
	}
	pjulian ju;
	ju = (pjulian)malloc(sizeof(JULIANDAY));
	utctojulianday(ut,ju);
	juliandaytobdt(ju, bt);
	free(ju);
}

void utctogpst(ptc ut, pgpst gt)
{
	pjulian ju;
	ju = (pjulian)malloc(sizeof(JULIANDAY));
	utctojulianday(ut, ju);
	juliandaytogpst(ju, gt);
	free(ju);
}

void utctojulianday(ptc ut, pjulian ju)
{
	int		m;
	int		y;
	double	dhour;

	dhour = ut->hour + ut->minute / (double)_HOUR_IN_MINUTE
		+ ut->second / (double)_HOUR_IN_SECOND;

	if (ut->month <= 2) {
		y = ut->year - 1;
		m = ut->month + 12;
	}
	else {
		y = ut->year;
		m = ut->month;
	}

	ju->daynum = (long)(365.25*y) + (long)(30.6001*(m + 1))
		+ ut ->day + (long)(dhour / 24 + 1720981.5);
	ju->secondnum = ((ut->hour + 12) % _DAY_IN_HOUR)*_HOUR_IN_SECOND
		+ ut->minute*_MINUTE_IN_SECOND + (long)ut->second;
	ju->secondfrac = ut->second-(long)ut->second;
}

void juliandaytoutc(pjulian ju, ptc ut)
{
	int a, b, c, d, e;
	double JD;
	JD = ju->daynum + (ju->secondnum + ju->secondfrac) / _DAY_IN_SECOND;

	a = static_cast<int>(JD + 0.5);
	b = a + 1537;
	c = static_cast<int>((b - 122.1) / 365.25);
	d = static_cast<int>(365.25*c);
	e = static_cast<int>((b - d) / 30.6001);
	
	double day = b - d-(long)(30.6001*e) + JD + 0.5 - a;
	ut->day = int(day);
	ut->month = e - 1 - 12 * (int)(e / 14);
	ut->year = c - 4715 - (int)((7 + ut->month) / 10);

	ut->hour = int((day - ut->day)*24.0);
	ut->minute = (int)(((day - ut->day)*24.0 - ut->hour)*60.0);
	ut->second = ju->secondnum + ju->secondfrac - (int((ju->secondnum + ju->secondfrac)/60))*60.0;

}

void gpsttojulianday(pgpst gt, pjulian ju)
{
	double JD;
	JD = gt->weeknum * 7 + (gt->secondnum + gt->secondfrac) / _DAY_IN_SECOND + 2444244.5;
	ju->daynum= long(JD);
	
	ju->secondnum = long(gt->secondnum + (gt->weeknum * 7 + 2444244.5 - ju->daynum)*_DAY_IN_SECOND);
	ju->secondfrac = gt->secondfrac;
}


void juliandaytogpst(pjulian ju, pgpst gt)
{
	double JD;
	JD = ju->daynum +( ju->secondnum + ju->secondfrac) / _DAY_IN_SECOND;
	gt->weeknum = int((JD - 2444244.5) / 7);

	gt->secondnum = long((JD - 2444244.5 - gt->weeknum * 7)*_DAY_IN_SECOND);
	gt->secondfrac = ju->secondfrac;

}

//�������ղ�ֵ
double deltjulianday(ptc u1, ptc u2)
{
	JULIANDAY j1, j2;
	utctojulianday(u1, &j1);
	utctojulianday(u2, &j2);
	double delt,d1,d2;
	d1 = j1.daynum + (j1.secondnum + j1.secondfrac) / _DAY_IN_SECOND;
	d2 = j2.daynum + (j2.secondnum + j2.secondfrac) / _DAY_IN_SECOND;
	//delt = (ju1->daynum - ju2->daynum)*_DAY_IN_SECOND + (ju1->secondnum - ju2->secondnum) + (ju1->secondfrac - ju2->secondfrac);
	delt = (d1 - d2)*_DAY_IN_SECOND;

	/*if (delt>302400)
		delt -= 604800;
	else if (delt<-302400)
		delt += 604800;
	else
		delt = delt;*/
	return delt;
}

//�����ձ任
void transjulian(pjulian ju1, double * dt, pjulian ju2)
{
	double JDold, JDnew;
	JDold = ju1->daynum + (ju1->secondnum + ju1->secondfrac) / _DAY_IN_SECOND;
	JDnew = JDold-*dt / _DAY_IN_SECOND;

	ju2->daynum = long(JDnew);
	ju2->secondnum = long((JDnew - long(JDnew))*_DAY_IN_SECOND);
	ju2->secondfrac = (JDnew - long(JDnew))*_DAY_IN_SECOND
		- long((JDnew - long(JDnew))*_DAY_IN_SECOND);
}






//���꺯��

//�ռ�ֱ������ϵ���������ϵ
void xyztoblh(pxyz px, pblh pb)
{
	double pi = 4.0*atan(1.0);
	double E2 = 2.0*flattening - flattening * flattening;
	double E4 = E2*E2;
	double ALFA = (px->x*px->x + px->y*px->y + (1.0 - E2)*px->z*px->z) / (a*a);
	double BATA = (px->x*px->x + px->y*px->y - (1.0 - E2)*px->z*px->z) / (a*a);
	double Q = 1.0 + 13.50*E4*(ALFA*ALFA - BATA*BATA) / pow(ALFA - E4, 3);
	double A1 = -Q + sqrt(Q*Q - 1.0);
	double AL = (1.0 / 3.0)*log(-A1);
	AL = -exp(AL);
	double A2 = AL + 1.0 / AL;
	double A3 = AL - 1.0 / AL;
	double T23 = (ALFA + E4 / 2.0) / 3.0 - (ALFA - E4)*A2 / 12.0;
	double T32 = sqrt(T23 *T23 + ((ALFA - E4)*A3)*((ALFA - E4)*A3) / 48.0);
	double T1 = -E2*BATA / (4.0*T32);
	double DK = T1 + sqrt(T23 + T32) - (1.0 - E2 / 2.0);
	double EK = (1.0 + DK) / (1.0 + DK - E2);
	/*cout << "T1��ֵ��" << T1 << endl;
	cout << "T32��ֵ��" << T32 << endl;
	cout << "T23��ֵ��" << T23 << endl;
	cout << "E2��ֵ��" << E2 << endl;
	cout << "ALFA��ֵ��" << ALFA << endl;
	cout << "AL��ֵ��" << AL << endl;
	cout << "DK��ֵ��" << DK << endl;
	cout << "EK��ֵ��" << EK << endl;*/
	pb->height = (DK / (1.0 + DK))*sqrt(px->x*px->x + px->y*px->y + (EK*px->z)*(EK*px->z));
	double P = sqrt(px->x*px->x + px->y*px->y);
	pb->latitude = atan(EK*px->z / P);
	double COSFL = px->x / P;
	double SINFL = px->y / P;
	pb->longitude = asin(SINFL);
	if (SINFL > 0.0&&COSFL<0.0) pb->longitude = pi - pb->longitude;
	if (SINFL < 0.0&&COSFL>0.0) pb->longitude = 2.0*pi + pb->longitude;
	if (SINFL <0.0&&COSFL<0.0) pb->longitude = pi - pb->longitude;

	/*double e2;//��һƫ���ʵ�ƽ��
	e2 = 2 * flattening - flattening*flattening;

	pb->longitude = atan(px->y / px->x);
	double W, N, N1 = 0, B, B1;
	B1 = atan(px->z / sqrt(px->x*px->x + px->y*px->y));
	while (1)
	{
		W = sqrt(1 - e2*sin(B1)*sin(B1));
		N1 = a / W;
		B = atan((px->z + N1*e2*sin(B1)) / sqrt(px->x*px->x + px->y*px->y));

		if (fabs(B - B1)<delta)
			break;
		else
			B1 = B;
	}

	pb->latitude = B;
	N = a / sqrt(1 - e2*sin(pb->latitude)*sin(pb->latitude));
	pb->height = sqrt(px->x*px->x + px->y*px->y) / cos(B) - N;*/
}

//�������ϵ���ռ�ֱ������ϵ
void blhtoxyz(pblh pb, pxyz px)
{
	double e2;//��һƫ���ʵ�ƽ��
	double N;//î��Ȧ�뾶
	e2 = 2 * flattening - flattening*flattening;
	N = a / sqrt(1 - e2*sin(pb->latitude)*sin(pb->latitude));

	px->x = (N + pb->height)*cos(pb->latitude)*cos(pb->longitude);
	px->y = (N + pb->height)*cos(pb->latitude)*sin(pb->longitude);
	px->z = (N*(1 - e2) + pb->height)*sin(pb->latitude);
}

//�ѿ�������ϵ��վ�Ŀռ�ֱ������ϵ
void xyztoenu(pxyz pxcenter, pxyz px, penu pe)
{
	double dx, dy, dz;
	dx = px->x - pxcenter->x;
	dy = px->y - pxcenter->y;
	dz = px->z - pxcenter->z;

	pblh pd;
	pd = (pblh)malloc(sizeof(BLH));

	xyztoblh(pxcenter, pd);

	pe->northing = -sin(pd->latitude)*cos(pd->longitude)*dx
		- sin(pd->latitude)*sin(pd->longitude)*dy
		+ cos(pd->latitude)*dz;
	pe->easting = -sin(pd->longitude)*dx
		+ cos(pd->longitude)*dy;
	pe->upping = cos(pd->latitude)*cos(pd->longitude)*dx
		+ cos(pd->latitude)*sin(pd->longitude)*dy
		+ sin(pd->latitude)*dz;
	free(pd);
}

//վ�Ŀռ�ֱ������ϵ���ѿ�������ϵ
 void enutoxyz(pxyz pxcenter, penu pe, pxyz px)
{
	pblh pd;
	pd = (pblh)malloc(sizeof(BLH));
	xyztoblh(pxcenter, pd);
	MatrixXd H(3, 3), DB(3, 1), DX(3, 1);
	DB(0, 0) = pe->northing;
	DB(1, 0) = pe->easting;
	DB(2, 0) = pe->upping;
	H(0, 0) = -sin(pd->latitude)*cos(pd->longitude);
	H(0, 1) = -sin(pd->latitude)*sin(pd->longitude);
	H(0, 2) = cos(pd->latitude);
	H(1, 0) = -sin(pd->longitude);
	H(1, 1) = cos(pd->longitude);
	H(1, 2) = 0;
	H(2, 0) = cos(pd->latitude)*cos(pd->longitude);
	H(2, 1) = cos(pd->latitude)*sin(pd->longitude);
	H(2, 2) = sin(pd->latitude);
	DX = (H.inverse())*DB;
	double dx, dy, dz;
	dx = DX(0,0 );
	dy = DX(1, 0);
	dz = DX(2, 0);
	px->x = pxcenter->x + dx;
	px->y = pxcenter->y + dy;
	px->z = pxcenter->z + dz;
	free(pd);
}

 //վ��ֱ������ϵ��վ�ļ�����ϵ
 void enutoenupolar(penu pe, penupolar pep)
 {
	 pep->range = sqrt(pe->northing*pe->northing + pe->easting*pe->easting + pe->upping*pe->upping);
	 pep->azimuth = atan(pe->easting / pe->northing);
	 //atan2�������صķ�ΧΪ(-PI,PI]��������ֵ������ʱ����1,2���ޣ�С����ʱ��3,4����
	 pep->azimuth = atan2(pe->easting / sqrt(pe->northing*pe->northing + pe->easting*pe->easting), pe->northing / sqrt(pe->northing*pe->northing + pe->easting*pe->easting));
	 if (pep->azimuth < 0.0)pep->azimuth += PI*2.0;
	 pep->elevation = atan(pe->upping / sqrt(pe->northing*pe->northing + pe->easting*pe->easting));
 }
 //�����Ǹ߶ȽǺͷ�λ��
 void sate_azi_ele(pxyz pxcenter, pxyz px, penupolar pep){
	 penu pe;
	 pe = (penu)malloc(sizeof(ENU));
	 xyztoenu(pxcenter,  px, pe);
	 enutoenupolar(pe, pep);
	 //cout << "���Ƿ�λ�ǣ�" << pep->azimuth*180.0/PI << endl;
	 //cout << "���Ǹ߶Ƚǣ�" << pep->elevation*180.0 / PI << endl;
	 //cout << "���Ǿ��룺" << pep->range<< endl;
	 free(pe);
 }
 //�󴩴̵���������Լ�ͶӰ��
 void ipp_pos(pblh pb, penupolar pep,pblh pb1,double* MF){
	 double psi;//�Ž�
	 psi = PI / 2.0 - pep->elevation - asin(ave_a / (ave_a + hion)*cos(pep->elevation));
	 pb1->latitude = asin(sin(pb->latitude)*cos(psi) + cos(pb->latitude)*sin(psi)*cos(pep->azimuth));
	 pb1->longitude = pb->longitude + atan(cos(pb->latitude)*sin(psi)*sin(pep->azimuth) / (cos(psi) - sin(pb->latitude)*sin(pb1->latitude)));
	 if (pb1->longitude < 0.0)pb1->longitude += PI*2.0;
	 if (pb1->longitude > PI*2.0)pb1->longitude -= PI*2.0;
	 pb1->height = hion;
	 *MF = 1.0 / sqrt(1.0 - (ave_a / (ave_a + hion)*cos(pep->elevation)*ave_a / (ave_a + hion)*cos(pep->elevation)));
	 //cout << "���̵�γ�ȣ�" << pb1->latitude*180.0 / PI << endl;
	 //cout << "���̵㾭�ȣ�" << pb1->longitude*180.0 / PI << endl;
	 //cout << "��б�ǣ�" << *MF << endl;
 }
 //�������ת���չ̵ش�����
 void g2m(pblh pb,double mjd,double pole_lat,double pole_lon)
 {
	 double mag_lat, mag_lon,sun_lon;
	 //�������ת���ش�����
	 mag_lat = asin(sin(pole_lat)*sin(pb->latitude) + cos(pole_lat)*cos(pb->latitude)*cos(pb->longitude - pole_lon));
	 //mag_lon = atan(cos(pole_lat)*cos(pb->latitude)*sin(pb->longitude - pole_lon) / (sin(mag_lat)*sin(pole_lat) - sin(pb->latitude)));
	 //if (fabs(fabs(mag_lon) - PI / 2.0) < 1e-8)cout <<"�ӽ��ٽ�ֵ��"<< mag_lon << endl;
	 mag_lon = atan2(cos(pb->latitude)*sin(pb->longitude - pole_lon) / cos(mag_lat), (sin(mag_lat)*sin(pole_lat) - sin(pb->latitude)) / (cos(mag_lat)*cos(pole_lat)));
	 pb->height = mag_lon;
	 //if (mag_lon < 0.0)mag_lon += PI*2.0;
	 cout.precision(10);
	 //cout << "���γ�ȣ�" << pb->latitude << endl;
	 //cout << "��ؾ��ȣ�" << pb->longitude << endl;
	 //cout << "�ش�γ�ȣ�" << mag_lat << endl;
	 //cout << "�شž��ȣ�" << mag_lon << endl;
	 //�ش�����ת���չ̵ش�����
	 sun_lon = PI*(1.0- 2.0 * (mjd - int(mjd)));//ƽ̫��������
	 //cout << "�����У�" << atan(sin(sun_lon - pole_lon) / sin(pole_lat) / cos(sun_lon - pole_lon)) << endl;
	 pb->latitude = mag_lat;
	 //double ang = atan2(sin(sun_lon - pole_lon) / sin(pole_lat), cos(sun_lon - pole_lon));
	 //if (ang < 0.0)ang += PI*2.0;
	 //pb->longitude = mag_lon - ang;
	 //pb->longitude = mag_lon - atan(sin(sun_lon - pole_lon) / sin(pole_lat)/cos(sun_lon - pole_lon));
	 pb->sun_lon = atan2(sin(sun_lon - pole_lon) / sin(pole_lat), cos(sun_lon - pole_lon));
	 pb->longitude = mag_lon - atan2(sin(sun_lon - pole_lon)/ sin(pole_lat), cos(sun_lon - pole_lon));
	 //pb->longitude = mag_lon - sun_lon;
	 if (pb->longitude < 0.0)pb->longitude += PI*2.0;
	 if (pb->longitude > PI*2.0)pb->longitude -= PI*2.0;
 }
 //��ionex�ļ���ֵ����������vtec
 double ionex_vtec(pio ionfile, UTC vtime, double lat, double lon){//lat��lon�Ǵ�ؾ�γ��
	 double p, q;//Ȩ
	 double v,v1, v2;//ǰ������ʱ���Ĳ�ֵ��С
	 //��ؾ�γ��תΪ���ľ�γ��
	 double lat_geo, lon_geo;
	 lon_geo = lon;
	 if (lon_geo > 180.0)lon_geo -= 360.0;
	 lat_geo = atan((1 - flattening)*(1 - flattening)*tan(lat*PI/180.0))*180.0/PI;
	 if (fabs(lon_geo) > 180.0 || fabs(lat_geo) > 87.5 || deltjulianday(&vtime, &ionfile->vtec_map[0].vtime)<0.0){
		 cout << "���̵㲻�ڼ��㷶Χ��" << endl;
		 return 0.0;
	 }
	 for (int i = 0; i < ionfile->vtec_map.size(); i++){
		 if (deltjulianday(&vtime, &ionfile->vtec_map[i].vtime) <= 0.0)//�ҵ���ֵʱ���������Ԫ
		 {
			 for (int j = 0; j < 71; j++)
			 {
				 if (lat_geo > 87.5 - j*2.5)//�ҵ���ֵ��γ�Ⱥ����γ�Ƚڵ�
				 {
					 for (int k = 0; k < 73; k++)
					 {
						 if (lon_geo < -180.0 + k *5.0)//�ҵ���ֵ�㾭�Ⱥ���ľ��Ƚڵ�
						 {
							 p = fabs((lat_geo - (87.5-j*2.5))/2.5);
							 q = fabs((lon_geo -(-180.0+(k-1)*5.0)) / 5.0);
							//�ռ��ֵ
							 //ǰһʱ���
							 v1 = (1 - p)*(1 - q)*ionfile->vtec_map[i-1].tec_values[j][k - 1] +
								  (1 - p)*q*ionfile->vtec_map[i-1].tec_values[j][k] +
								   p*(1 - q)*ionfile->vtec_map[i-1].tec_values[j-1][k - 1] +
								   p*q*ionfile->vtec_map[i-1].tec_values[j-1][k];
							 //��һʱ���
							 v2 = (1 - p)*(1 - q)*ionfile->vtec_map[i].tec_values[j][k - 1] +
								 (1 - p)*q*ionfile->vtec_map[i].tec_values[j][k] +
								 p*(1 - q)*ionfile->vtec_map[i].tec_values[j - 1][k - 1] +
								 p*q*ionfile->vtec_map[i].tec_values[j - 1][k];
							/* if (vtime.hour == 1 && vtime.minute == 59 && vtime.second > 0.0){
								 cout << "1:" << "i:" << i << "j:" << j << "k:" << k << "lat:" << lat_geo << "lon:" << lon_geo<< endl;
								 cout << ionfile->vtec_map[i - 1].tec_values[j][k - 1] <<" "<< ionfile->vtec_map[i - 1].tec_values[j][k]
									 << " " << ionfile->vtec_map[i - 1].tec_values[j - 1][k - 1] << " " << ionfile->vtec_map[i - 1].tec_values[j - 1][k] << " " << v1 << " " << v2<<endl;
								 cout << ionfile->vtec_map[i ].tec_values[j][k - 1] << " " << ionfile->vtec_map[i ].tec_values[j][k]
									 << " " << ionfile->vtec_map[i ].tec_values[j - 1][k - 1] << " " << ionfile->vtec_map[i ].tec_values[j - 1][k] << " " << v1 << " " << v2 << endl;
							 }
							 if (vtime.hour == 2 && vtime.minute == 0 && vtime.second <30.0){
								 cout << "2:" << "i:" << i << "j:" << j << "k:" << k<< "lat:" << lat_geo << "lon:" << lon_geo << endl;
								 cout << ionfile->vtec_map[i - 1].tec_values[j][k - 1] << " " << ionfile->vtec_map[i - 1].tec_values[j][k]
									 << " " << ionfile->vtec_map[i - 1].tec_values[j - 1][k - 1] << " " << ionfile->vtec_map[i - 1].tec_values[j - 1][k] << " " << v1 << " " << v2 << endl;
								 cout << ionfile->vtec_map[i].tec_values[j][k - 1] << " " << ionfile->vtec_map[i].tec_values[j][k]
									 << " " << ionfile->vtec_map[i].tec_values[j - 1][k - 1] << " " << ionfile->vtec_map[i].tec_values[j - 1][k] << " " << v1 << " " << v2 << endl;
							 }*/
							 //ʱ���ֵ
							 v = (fabs(deltjulianday(&vtime, &ionfile->vtec_map[i - 1].vtime))*v2 +
								 fabs(deltjulianday(&vtime, &ionfile->vtec_map[i].vtime))*v1) / fabs(deltjulianday(&ionfile->vtec_map[i].vtime, &ionfile->vtec_map[i-1].vtime));
							 return v;
						 }
					 }
				 }
			 }
		 }
	 }
	 cout << "û���ҵ���ֵʱ���" << endl;
	 return 0.0;
 }
 //����黯���õ¶���ʽ��г��
 void geo_legendre(double lat, double**pg,int n)
 {
	 double sin_lat = sin(lat);
	 double tmp,tmp1,tmp2;
	 pg[0][0] = sqrt(3.0 * (1.0 - sin_lat*sin_lat));
	 //����������г��
	for (int i = 1; i < n; i++){
		tmp = sqrt((2.0 * i + 3.0) / (2.0 * i + 2.0)*(1.0 - sin_lat*sin_lat));
		pg[i][i] = tmp*pg[i - 1][i - 1];
	}
	//�������һ���������г��
	for (int j = 1; j < n - 1; j++){
		for (int i = j + 1; i < n; i++){
			tmp1 =sqrt((2.0 * i + 3.0)*(2.0 * i + 1.0)/(i+j+2)/(i-j))*sin_lat;
			tmp2 =sqrt((2.0 * i + 3.0)*(i + j + 1)*(i - j - 1) / (2.0 * i - 1.0) / (i + j + 2.0) / (i - j));
			pg[i][j] = tmp1*pg[i-1][j] - tmp2*pg[i-2][j];
		}
	}
	//������г��ĵ�һ��
	int j = 0;
	for (int i = 1; i < n; i++)
	{
		tmp1 = sqrt(1.0*(2 * i + 3)*(2 * i + 1) / (i + j + 2) / (i - j))*sin_lat;
		tmp2 = sqrt(1.0*(2 * i + 3)*(i + j + 1)*(i - j - 1) / (2 * i - 1) / (i + j + 2) / (i - j));
		if (i == 1){ pg[i][j] = tmp1*pg[i - 1][j]; }
		else{ pg[i][j] = tmp1*pg[i - 1][j] - tmp2*pg[i - 2][j];}
	}
 
 }
 //����黯���õ¶���ʽ��г��
 void zone_legendre(double lat, double *pz, int n){
	 double sin_lat = sin(lat);
	 pz[0] = 1.0;
	 pz[1] = sqrt(3.0)*sin_lat;
	 double tmp1, tmp2, tmp3;
	 for (int i = 2; i < n + 1; i++)
	 {
		 tmp1 = (2.0 - 1.0 / i)*sin_lat;
		 tmp2 = sqrt((2.0 * i - 1.0) / (2.0* i - 3.0))*(1.0 - 1.0 / i);
		 tmp3 = sqrt((2.0 * i + 1.0) / (2.0 *i - 1.0));
		 pz[i] = tmp3*(tmp1*(pz[i - 1]) - tmp2*pz[i - 2]);
	 }
 }
 //����黯���õº�����ֵ
 void legendre(double lat, double lon, double* legend,int n){//��(1+2L+1)(O+1)/2=(O+1)^2��
	 double**pg = new double*[n];
	 for (int i = 0; i < n; i++){
		 pg[i] = new double[n];
	 }
	 double* pz = new double[n + 1];
	 for (int i = 0; i < n; i++){//��ʼ��
		 pz[i] = 0.0;
		 for (int j = 0; j < n; j++){
			 pg[i][j] =0.0;
		 }
	 }
	 pz[n] = 0.0;
	 //double pg[O][O] = { 0.0 };
	// double pz[O+1] = { 0.0 };
	 geo_legendre(lat, pg,n);
	 zone_legendre(lat, pz,n);
	 cout.setf(ios_base::fixed, ios_base::floatfield);
	 cout.precision(10);
	 for (int i = 0; i < n + 1; i++)
	 {
		 legend[i*i] = pz[i];
		 //cout << "pz[" << i << "]:" << pz[i] << endl;
		 //if (i == O){ break; }
		 for (int j = 0; j < i; j++)
		 {
			legend[i*i+ j*2+1] = pg[i-1][j] * cos((j+1)*lon);
			legend[i*i+ j*2+2] = pg[i-1][j] * sin((j+1)*lon);
			if (legend[i*i + j * 2 + 1] == 0.0 || legend[i*i + j * 2 + 2] == 0.0){ cout << "                                               ���õ�ֵcʧ��" << endl;}
			//cout << "pg[" << i-1 << "]" << "[" << j << "]:" << pg[i-1][j] << endl;
		 }
	 }
	 /*cout << "γ�ȣ�:" <<lat << endl;
	 for (int i = 0; i < O + 1; i++){
		 cout << "legend[" << i << "]:" << legend[i] << endl;
	 }*/
	 for (int i = 0; i < n; i++){
		 delete[] pg[i];
	 }
	 delete[] pg;
	 delete pz;
 }
 //����г�����ļ���ֵ����������vtec
 double sh_vtec(psh sh_file, UTC vtime, double lat, double lon,int n){//lat��lon���չ̵ش�����ϵ������
	 double *legend = new double[(n + 1)*(n + 1)];
	 double mjd;
	 double  dt1, dt2;//dt1������ʱ���ǰһ��Ԫ�Ĳ�ֵ��ֵ��dt2������ʱ��ͺ�һ��Ԫ�Ĳ�ֵ��ֵ
	 JULIANDAY jul;
	 one_coeff coe1,coe2;
	 BLH ipp;
	 //cout << "�۲�ʱ�䣺" << vtime.hour << ":" << vtime.minute << ":" << vtime.second << endl;
	 ipp.latitude = lat; ipp.longitude = lon; ipp.height = hion;
	 double vtec = 0.0;
	 utctojulianday(&vtime, &jul);//�Ѽ��㴩�̵�ʱ��utcת��������
	 mjd = jul.daynum + (jul.secondnum + jul.secondfrac) / 86400.0 - 2400000.5;//������ת��Լ��������
	 for (int i = 0; i < sh_file->coeff.size(); i++)
	 {
		 if (deltjulianday(&vtime, &sh_file->coeff[i].jul) <= 0.0 || i == (sh_file->coeff.size()-1))//�ҵ���ֵʱ���������Ԫ�������һ��
		 {//��ʱû��ʹ��ʱ���ֵ
			 if (i == 0){ coe1 = sh_file->coeff[i]; coe2 = sh_file->coeff[i]; }//��һ����Ԫ
			 else if (i == (sh_file->coeff.size() - 1) && deltjulianday(&vtime, &sh_file->coeff[i].jul) > 0.0){ coe1 = sh_file->coeff[i]; coe2 = sh_file->coeff[i]; }//���һ����Ԫ����������ʱ�䲻�Ǻ�һ������
			 else{ coe1 = sh_file->coeff[i - 1]; coe2 = sh_file->coeff[i-1]; }//��г����ϵ������Ч��Χһ��Ϊ���1����2��Сʱ������ȡ�۲�ʱ��ǰ�����Ԫ
			 //dt1 = deltjulianday(&vtime, &coe1.jul) / deltjulianday(&coe2.jul, &coe1.jul); dt2 = 1.0 - dt1;
			 dt1 = 1.0; dt2 = 1.0 - dt1;
			 //cout << "shʱ�䣺" << coe1.jul.hour << ":" << coe1.jul.minute << ":" << coe1.jul.second << endl;
			 //cout << "2ʱ�䣺" << coe2.jul.hour << ":" << coe2.jul.minute << ":" << coe2.jul.second << endl;
			 cout.setf(ios::fixed);
			 cout.precision(10);
			 //cout << "���̵�γ�ȣ�" << ipp->latitude<< endl;
			 //cout << "���̵㾭�ȣ�" << ipp->longitude<< endl;
			// cout << "Լ�������գ�" << mjd << endl;
			 g2m(&ipp, mjd, coe1.pole_lat*PI / 180.0, coe1.pole_lon*PI / 180.0);//���̵��ɴ������ת���չ̵ش�����
			 //cout << "�չ̵ش�γ�ȣ�" << ipp->latitude << endl;
			 //cout << "�չ̵شž��ȣ�" << ipp->longitude << endl;
			 //cout << "���̵�γ������ֵ��" <<sin(ipp->latitude) << endl;
			 legendre(ipp.latitude, ipp.longitude, legend,n);//�ɴ��̵�������õ¶���ʽ��ֵ
			// cout << "���õ�ֵ" << endl;
			//for (int k = 0; k < 256; k++){
				 //cout <<"legend["<<k<<"]:"<< legend[k] << endl;
			 //}
			// cout << "���̵��չ̵ش����꾭γ�ȣ�" << ipp.latitude << ":" << ipp.longitude<< endl;
			 //cout << "ϵ��ʱ�䣺" << coe1.jul.hour << ":" << coe1.jul.minute << ":" << coe1.jul.second << endl;
			 for (int j = 0; j <coe1.one_coe.size(); j++)
			 {
				 //cout << coe1.one_coe[j].degree << "��" << coe1.one_coe[j].order << "��"<<"ϵ����" << coe1.one_coe[j].tec << endl;
				 //cout << "���õ�ϵ����" << legend[j] << endl;
				 vtec += (coe1.one_coe[j].tec*dt2 + coe2.one_coe[j].tec*dt1) * legend[j];
				 //vtec += coe1.one_coe[j].tec* legend[j];
			 }
			 return vtec;
			 delete legend;
		 }
	 }
	 cout << "û���ҵ���Ч����г����ϵ��" << endl; 
	 delete legend;
	 return 0.0;
 }
 //����ת��
 void Transposition(double**L, int n)
 {
	 for (int i = 0; i < n; i++)
	 {
		 for (int j = 0; j < i; j++)
			 swap(L[i][j], L[j][i]);
	 }
 }
 //����˷�
 void Multi(double**A, double**B, int n)//AXB->B
 {
	 double **C = new double*[n];
	 for (int i = 0; i < n; i++)
		 C[i] = new double[n];
	 for (int i = 0; i < n; i++)
	 {
		 for (int j = 0; j < n; j++)
		 {
			 C[i][j] = 0;
			 for (int k = 0; k < n; k++)
				 C[i][j] += A[i][k] * B[k][j];
		 }
	 }
	 for (int i = 0; i < n; i++)
	 {
		 for (int j = 0; j < n; j++)
			 B[i][j] = C[i][j];
	 }
	 for (int i = 0; i < n; i++)
	 {
		 delete[] C[i];
		 C[i] = NULL;
	 }
	 delete C;
	 C = NULL;
 }
 //Cholesky�ֽⷽ��
 void chol_rf(double**A, double**L, double**d, int n)
 {
	 double s1 = 0.0;
	 double **g = new double*[n]; //������  
	 for (int i = 0; i < n; i++)
		 g[i] = new double[n]; //������  
	 for (int i = 0; i < n; i++){//��ʼ��
		 for (int j = 0; j < n; j++){
			 g[i][j] = 0.0;
		 }
	 }
	 d[0][0] = A[0][0];
	 for (int i = 1; i < n; i++){
		 for (int j = 0; j < i; j++){
			 s1 = 0.0;
			 for (int k = 0; k < j; k++){
				 s1 += g[i][k] * L[j][k];
			 }
			 g[i][j] = A[i][j] - s1;
		 }
		 for (int j = 0; j < i; j++){
			 L[i][j] = g[i][j] / d[j][j];
		 }
		 s1 = 0.0;
		 for (int k = 0; k < i; k++){
			 s1 += g[i][k] * L[i][k];
		 }
		 d[i][i] = A[i][i] - s1;
	 }
	 for (int i = 0; i < n; i++){
		 L[i][i] = 1.0;
	 }
	 for (int i = 0; i < n; i++)
	 {
		 delete[] g[i];
	 }
	 delete[] g;
 }
 /*
 void chol_rf(double**A, double**L, double**d, int n)
 {
	 double s1 = 0.0, s2 = 0.0;;
	 for (int j = 0; j < n; j++)
	 {
		 s1 = 0.0;
		 for (int k = 0; k < j; k++)
		 {
			 s1 += L[j][k] * L[j][k] * d[k][k];
		 }
		 d[j][j] = A[j][j] - s1;
		 for (int i = j + 1; i < n; i++)
		 {
			 s2 = 0.0;
			 for (int k = 0; k < j; k++)
			 {
				 s2 += L[i][k] * L[j][k] * d[k][k];
			 }
			 L[i][j] = (A[i][j] - s2) / d[j][j];
		 }
		 L[j][j] = 1.0;
	 }
 }*/
 //Cholesky�ֽⷽ���ⷽ����
 void chol_eq(double**A, double* b, double* x,int n)
 {
	 //cout << "Cholesky�ֽⷽ���ⷽ����" << endl;
	 double *y = new double[n];
	 double **L = new double*[n]; //������  
	 double **d = new double*[n]; //������
	 for (int i = 0; i < n; i++)
	{	L[i] = new double[n]; //������ 
		d[i] = new double[n]; //������ 
	}
	 for (int i = 0; i < n; i++){//��ʼ��
		 y[i] = 0.0;
		 for (int j = 0; j < n; j++){
			 L[i][j] = 0.0;
			 d[i][j] = 0.0;
		 }
	 }
	 chol_rf(A, L, d,n);//Cholesky�ֽ�
	 y[0] = b[0];
	 double s1;
	 for (int i = 1; i < n; i++){
		 s1 = 0.0;
		 for (int k = 0; k < i; k++){
			 s1 += L[i][k] * y[k];
		 }
		 y[i] = b[i] - s1;
	 }
	 x[n - 1] = y[n - 1] / d[n - 1][n - 1];
	 for (int i = n - 2; i >= 0; i--){
		 s1 = 0.0;
		 for (int k = i+1; k < n; k++){
			 s1 += L[k][i] * x[k];
		 }
		 x[i] = y[i] / d[i][i] - s1;
	 }
	 for (int i = 0; i < n; i++)
	 {
		 delete[] L[i];
		 delete[] d[i];
	 }
	 delete y;
	 delete[] L;
	 delete[] d;
 }
 /*
 void chol_eq(double**A, double* b, double* x,int n)
 {
	 //cout << "Cholesky�ֽⷽ���ⷽ����" << endl;
	 double **L = new double*[n]; //������  
	 double **d = new double*[n]; //������
	 for (int i = 0; i < n; i++)
	 {
		 L[i] = new double[n]; //������ 
		 d[i] = new double[n]; //������ 
	 }
	 for (int i = 0; i < n; i++){//��ʼ��
		 for (int j = 0; j < n; j++){
			 L[i][j] = 0.0;
			 d[i][j] = 0.0;
		 }
	 }
	 chol_rf(A, L, d,n);//Cholesky�ֽ�
	 
 for (int k = 0; k < n; k++)
 {
	 for (int i = 0; i < k; i++)
		 b[k] -= b[i] * L[k][i];
	 b[k] /= L[k][k];
 }
 
 Transposition(L, n);//L^T
 Multi(d, L, n);//D X L^T
 for (int k = n - 1; k >= 0; k--)
 {
	 for (int i = k + 1; i < n; i++)
		 b[k] -= b[i] * L[k][i];
	 b[k] /= L[k][k];
 }
 for (int k = 0; k < n; k++)
 {
	 x[k] = b[k];
 }
 for (int i = 0; i < n; i++)
 {
	 delete[] L[i];
	 delete[] d[i];
 }
 delete[] L;
 delete[] d;
 }
 */
 //Ѱ��chazh����
 bool find_sp3_ephem(string prn, ptc ut, psp3 sp3file, s_sp3_ephe sp3[])//n���������ղ�ֵ��������sp3�����صĲ�ֵ�ڵ�ֵ
 {
	 int num;//��ֵ��ǰ��Ľڵ���
	 if (NN % 2 == 0)num = NN / 2;//��ֵ����Ϊż��
	 else num = (NN + 1) / 2;
	// JULIANDAY jld1, jld2;//jlld1�����źŷ���ʱ�䣬jlld2����������ʱ��
	 /*cout << "�����źŷ���ʱ�䣺" << ut->year<< ":" << ut->month << ":" << ut->day;
	 cout << ":" << ut->hour << ":" << ut->minute << ":" << ut->second << endl; */
	 int kk = 32;
	 for (int i = 0; i < kk; i++)
	 {
		 if (prn.substr(0, 3) == sp3file->all_ephem[i].prn.substr(0, 3))
		 {
			 //cout << "�ҵ���Ӧ�����Ǻ�" << endl;
			 // satn = true;
			 int w = sp3file->all_ephem[i].sate_ephem.size();
			 for (int j = 0; j < w; j++)
			 {
				 double delta;
				 delta = deltjulianday(ut, &sp3file->all_ephem[i].sate_ephem[j].utime_n);
				 if (delta > 0.0) continue;//�ҵ���ֵ��ǰ��Ĳ�ָ��ʱ��
				 else{
					 //cout << "�ҵ���Ԫ" << endl;
					 if (j == 0) return false;//��ָ�������в�ֵ�ڵ�֮�⣨֮ǰ��
					 else if (j  < num){//j�պ��ǲ�ֵʱ���ǰ��Ŀ��ò�ֵ���������ӵ�һ����ֵ�㿪ʼȡǰn+1���ڵ�������ֵ
						 for (int k = 0; k < NN + 1; k++)
						 {
							 sp3[k] = sp3file->all_ephem[i].sate_ephem[k];
							// cout << "ǰ" << "x:" << sp3[k].x << " y:" << sp3[k].y << " z:" << sp3[k].z << endl;
						 }
					 }
					 else if (w - j < NN + 1 - num){//��ֵ��������棬����ڵ������㹻
						 for (int k = 0; k < NN + 1; k++)
						 {
							 sp3[k] = sp3file->all_ephem[i].sate_ephem[w-1 - NN + k];
							 //cout << "��" << "x:" << sp3[k].x << " y:" << sp3[k].y << " z:" << sp3[k].z << endl;
						 }
					 }
					 else{
						 for (int k = 0; k < NN + 1; k++)
						 {
							 sp3[k] = sp3file->all_ephem[i].sate_ephem[j - num + k];//j�ǵ�num+1��
							// cout << "��" << "x:" << sp3[k].x << " y:" << sp3[k].y << " z:" << sp3[k].z << endl;
						 }
					 }
					 return true;
				 }
			 }
		 }
	 }
	 // cout
	 return false;
 }
 
 //��������������������
 bool cal_sp3_sate_coor(string prn, ptc ut, psp3 sp3file, pxyz coor){
	 coor->x = 0.0;
	 coor->y = 0.0;
	 coor->z = 0.0;
	 s_sp3_ephe sp3[NN + 1];
	 double s = 1.0;
	 if(find_sp3_ephem(prn, ut, sp3file, sp3))
	 {
		 //cout <<"�ҵ�����"<< endl;
		 for (int i = 0; i < NN + 1; i++)
		 {//�������շ������������
			 s = 1.0;
			 for (int j = 0; j < NN + 1; j++)
			 {
				 if (i == j)continue;
				 //cout << "x-xj:" << deltjulianday(ut, &sp3[j].utime_n) << "xi-xj:" << deltjulianday(&sp3[i].utime_n, &sp3[j].utime_n) << endl;
				 s = s*deltjulianday(ut, &sp3[j].utime_n) / deltjulianday(&sp3[i].utime_n, &sp3[j].utime_n);
			 }
			 coor->x += s*sp3[i].x;
			 coor->y += s*sp3[i].y;
			 coor->z += s*sp3[i].z;
		 }
		 return true;
	 }
	 else{ 
		 return false; 
	 }

 }
 //��dcb�ļ����ҵ���Ӧ��dcbֵ
 bool find_dcb(bool b,string prn, string station, pdcb dcbfile, double* dcb_val){//b==0��������dcb��b!=0�����վdcb
	 //cout << "dcb�ļ���С��" << dcbfile->dcb_val.size() << endl;
	 for (int i = 0; i < dcbfile->dcb_val.size(); i++)
	 {
		 //cout << "����prn:" << dcbfile->dcb_val[i].prn << "  ��վ:" << dcbfile->dcb_val[i].station << "  obs1:" << dcbfile->dcb_val[i].obs1 << "  obs2:" << dcbfile->dcb_val[i].obs2 << endl;
		 //cout << "prn:" << prn << "  station:" << station <<endl;
		 if (((b==0&&dcbfile->dcb_val[i].prn == prn)||(b!=0&&dcbfile->dcb_val[i].station == station))&&dcbfile->dcb_val[i].obs1 == "C1C"&&dcbfile->dcb_val[i].obs2 == "C2W")
		 {
			 *dcb_val = dcbfile->dcb_val[i].bias*1e-9*c;//���ؾ���ֵ
			 return true;
		 }
	 }
	 return false;
 }
 //����c_a0
 double c_a0(double**ctable, UTC ut, double lat, double lon){
	 double legend[36] = { 0.0 };
	 double beta[17] = { 0.0 };//17��ͨ������������ϵ��
	 double b[17] = { 0.0 };//17���ʹ��̵��йص����õ�ϵ��ֵ
	 double a0 = 0.0;
	 JULIANDAY jul;
	 double mjd;
	 legendre(lat, lon, legend, 5);
	 //mjd ��Ӧ����Լ�������յ���������ʱ��
	 ut.minute = 0; ut.second = 0.0;
	 if (ut.hour % 2 == 0){ ut.hour += 1; }
	 utctojulianday(&ut, &jul);//�Ѽ��㴩�̵�ʱ��utcת��������
	 mjd = jul.daynum + (jul.secondnum + jul.secondfrac) / 86400.0 - 2400000.5; 
	 double per[12] = { 1.0, 0.5, 0.33, 14.6, 27.0, 121.6, 182.51, 365.25, 4028.71, 2014.35, 1342.9, 1007.18 };
	 for (int j = 0; j < 17; j++){
		 beta[j] = ctable[0][j];
		 for (int i = 0; i < 12; i++){
			 beta[j] += ctable[2 * i + 1][j] * cos(2 * PI / per[i] * mjd)+ctable[2 * i + 2][j] * sin(2 * PI / per[i] * mjd);
		 }
	 }
	 for (int j = 0; j < 12; j++){
		 b[j] = legend[9 + j];//3�����4����ǰ5��
	 }
	 for (int j = 12; j < 17; j++){
		 b[j] = legend[13 + j];//5����ǰ5��
	 }
	 for (int j = 0; j < 17; j++){
		 a0 += beta[j] * b[j];
	 }
	 return a0;
 }
 //������
 void putresult(pvt pt)
 {
	 ofstream outfile("E:\\GNSS\\���\\sh_vtec3.txt", ios::out);
	 if (!outfile)
	 {
		 cerr << "�ļ�����ʧ��" << endl;
		 abort();
	 }
	 cout << "��ʼд�ļ�" << endl;
	 for (unsigned int i = 0; i < pt->allresult.size(); i++)
	 {
		 outfile.setf(ios_base::fixed, ios_base::floatfield);
		 outfile.precision(8);
		 outfile <<setw(4) << pt->allresult[i].rtime.year << setw(3) << pt->allresult[i].rtime.month << setw(3) << pt->allresult[i].rtime.day << 
			 setw(3) << pt->allresult[i].rtime.hour << setw(3) << pt->allresult[i].rtime.minute << setw(15) << pt->allresult[i].rtime.second<<endl;
		 outfile << " �в�                                                                                      ������" << pt->allresult[i].residul.size()<< endl;
		 for (unsigned int j = 0; j < pt->allresult[i].residul.size(); j++){
			 outfile << std::right << setw(12) << pt->allresult[i].residul[j];
			 if ((j + 1) % 10 == 0)outfile << endl;
		 }
		 outfile << endl;
		 outfile << "COEFFICIENTS" << endl;
		 outfile << "DEGREE   ORDER        VALUE (TECU)        VALUE (TECU)" << endl;
		 for (unsigned int k = 0; k < O + 1; k++){
			 outfile << std::right << setw(4) << k << std::right << setw(9) << 0 << std::right << setw(20) << pt->allresult[i].coe[k*k] << std::right << setw(20) << pt->allresult[i].coe1[k*k] << endl;
			 for (int l = 1; l <= k; l++){
				 outfile << std::right << setw(4) << k << std::right << setw(9) << l << std::right << setw(20) << pt->allresult[i].coe[k*k + l * 2 - 1] << std::right << setw(20) << pt->allresult[i].coe1[k*k + l * 2 - 1] << endl;
				 outfile << std::right << setw(4) << k << std::right << setw(9) << -l << std::right << setw(20) << pt->allresult[i].coe[k*k + l * 2] << std::right << setw(20) << pt->allresult[i].coe1[k*k + l * 2] << endl;
			 }
		 }
	 }
	 outfile.close();
	 cout << "�ļ�д�����" << endl;
 }
 //�ɹ۲�ֵ����vtec
 double obs_vtec(pobs obsfile, UTC obs_t,int n, string prn, psp3 sp3file, pdcb dcbfile, pblh ipp)
 {
	 UTC day_first;//һ���0ʱ
	 day_first = obs_t;
	 day_first.hour = 0; day_first.minute = 0; day_first.second = 0.0;
	 int ns=0;//һ��֮�еĵ�n���۲���Ԫ,ns��Ҫ����Ĺ۲���Ԫ֮ǰ����ʼʱ��
	 int flag;//�ж��ϸ���Ԫ�Ƿ���Ȼ�۲⵽Ҫ���������
	 vector<int>prn_pos;//��i����Ԫ����prn��ŵ�λ��
	 int p;//λ��
	 int num = 0, qq = 0;//num:���̵�����,qq:�����۲�һ�����ǵĹ۲���Ԫ��,cycle_slip:��������
	 XYZ sta_coor;//��վwgs84����
	 BLH sta_coor_blh;//��վ�ʹ��̵�Ĵ������
	 pxyz sate_coor = new XYZ;//��������
	 penupolar pep = new ENUPOLAR;//���Ƿ�λ�Ǻ͸߶Ƚ�
	 double* mf = new double;//ͶӰ����
	 double dt2 = 0.0;//dt1�������źŴ���ʱ��,dt2��ͬһ������������Ԫ��ʱ���
	 double w = 1.0;//Ȩ
	 double mf0 = 1.0;//��һ����mf
	 double dcb_sat = 0.0, dcb_sta = 0.0;//��վ������dcb
	 //����̽��
	 double n_mw1 = 0.0, n_mw2 = 0.0;//MW���
	 double n0 = 0.0, n1 = 0.0;//��ֵ
	 double sigma0 = 0.0, sigma1 = 0.0;//��׼��
	 int cycle_slip = 0;//��������
	 double vtec1 = 0.0, vtec0 = 0.0;//��һ������õ���vtec0�͵�ǰ�����vtec1
	 if (!find_dcb(0, prn, obsfile->obsheaddata.station, dcbfile, &dcb_sat) || !find_dcb(1, prn, obsfile->obsheaddata.station, dcbfile, &dcb_sta))
	 {
		 cout << "û���ҵ�����" << prn << "���߲�վdcb" << endl;
		 delete sate_coor;
		 delete pep;
		 delete mf;
		 return 0.0;
	 }
	 sta_coor = obsfile->obsheaddata.approx_coordinate;
	 xyztoblh(&sta_coor, &sta_coor_blh);
	 //������ʼ�۲�ʱ��
	 for (int i = n; i >= 0; i--)
	 {
		 flag = 0;
		 for (int j = 0; j < obsfile->obsdata[i].prn_list.size(); j++)
		 {
			 if (prn == obsfile->obsdata[i].prn_list[j])
			 {
				 flag = 1;
				 prn_pos.push_back(j);
				 break;
			 }
		 }
		 if (flag == 0)//������Ҫ��������ǵĹ۲����ݾͽ�����
		 {
			 ns = i+1;
			 break;
		 }
		// if (i == 0)ns = 0;//��ʼʱ����һ��Ŀ�ʼ,��ʼ��ns=0����
	 }
	 //cout << "n:" << n<< endl;
	 //cout << "ns:" << ns << endl;
	 //cout << "λ�ø���:" << prn_pos.size() << endl;
	 //һ������������Ԫ������λƽ��α��ķ�����vtec
	 //int c = 0;
	 for (int j = ns; j <= n; j++)
	 {
		 //cout << "j:" << j << endl;
		 //cout << "qq:" << qq<< endl;
		 p = n - j;
		 if (!cal_sp3_sate_coor(prn, &obsfile->obsdata[j].utime_o, sp3file, sate_coor))//����������
		 {
			 //cout << "�����������ʧ��" << endl;
			 continue;
		 }//�����������ʧ��
		 sate_azi_ele(&sta_coor, sate_coor, pep);//�����Ƿ�λ�Ǻ͸߶Ƚ�
		 if (pep->elevation*180.0 / PI < 15.0){ continue; }//���Ǹ߶Ƚ�Ҫ����15��
		 ipp_pos(&sta_coor_blh, pep, ipp, mf);//�󴩴̵�γ�Ⱦ���,��λ����
		 if (j != ns && num>0)
		 {//��һ��������dt2
			 dt2 = fabs(deltjulianday(&obsfile->obsdata[j].utime_o, &obsfile->obsdata[j - 1].utime_o));
		 }
		 if (qq == 0 || num == 0 || dt2 > 300.0)//qq=0��ζ����һ�������������ߵ�һ����������Ԫʱ��������300s��֤������������������ҪѡΪʱ��������vtec
		 {
			/* cout << "1416yuejie" << endl;
			 cout << "p:" << p << endl;
			 cout << "prn_pos[p]:" << prn_pos[p] << endl; 
			 cout << "obsfile->obsdata[j].one_obs_data:" << obsfile->obsdata[j].one_obs_data.size() << endl;
			 cout << "obsfile->obsdata[j].prn_list:" << obsfile->obsdata[j].prn_list.size() << endl;
			 cout << "obsfile->obsdata[j].one_obs_data.prn_list[prn_pos[p]]:" << obsfile->obsdata[j].prn_list[prn_pos[p]] << endl;
			 cout << "obsfile->obsdata[j].one_obs_data.prn_list[prn_pos[p]]:" << obsfile->obsdata[j].one_obs_data[prn_pos[p]].p1 << endl;
			 vtec1 = 9.52437*(obsfile->obsdata[j].one_obs_data[prn_pos[p]].p2 - obsfile->obsdata[j].one_obs_data[prn_pos[p]].p1 + dcb_sta + dcb_sat) / (*mf);*/
			 w = 1.0;
			 //cout << "vtec1:" << vtec1 << endl;
			 n_mw1 = obsfile->obsdata[j].one_obs_data[prn_pos[p]].l1 - obsfile->obsdata[j].one_obs_data[prn_pos[p]].l2 - (f1 - f2) / (f1 + f2)*(obsfile->obsdata[j].one_obs_data[prn_pos[p]].p1*f1 / c + obsfile->obsdata[j].one_obs_data[prn_pos[p]].p2*f2 / c);
			 qq = 1;
			 n1 = n_mw1;
			 sigma1 = 0.0;
		 }
		 else{
			 n0 = n1;
			 sigma0 = sigma1;
			 //cout << "1426yuejie" << endl;
			 n_mw1 = obsfile->obsdata[j].one_obs_data[prn_pos[p]].l1 - obsfile->obsdata[j].one_obs_data[prn_pos[p]].l2 - (f1 - f2) / (f1 + f2)*(obsfile->obsdata[j].one_obs_data[prn_pos[p]].p1*f1 / c + obsfile->obsdata[j].one_obs_data[prn_pos[p]].p2*f2 / c);
			 //cout << "1427yuejie" << endl;
			 if (j<n){ n_mw2 = obsfile->obsdata[j + 1].one_obs_data[prn_pos[p - 1]].l1 - obsfile->obsdata[j + 1].one_obs_data[prn_pos[p - 1]].l2 - (f1 - f2) / (f1 + f2)*(obsfile->obsdata[j + 1].one_obs_data[prn_pos[p - 1]].p1*f1 / c + obsfile->obsdata[j + 1].one_obs_data[prn_pos[p - 1]].p2*f2 / c); }
			 n1 = n0 + 1.0 / qq*(n_mw1 - n0);
			 sigma1 = sqrt(sigma0*sigma0 + 1.0 / qq*((n_mw1 - n0)*(n_mw1 - n0) - sigma0*sigma0));
			 if (qq>2 && fabs(n_mw1 - n0) >= 4.0*sigma0&&fabs(n_mw2 - n_mw1) < 1.0){//�������������һ����Ԫn_mw2 == n_mw1����Ϊû�����¸�ֵ��n_mw2
				 cycle_slip++;
				 qq = 0;
				 continue;
			 }
			 w = 1.0 / qq; 
			 //cout << "1437yuejie" << endl;
			 vtec1 = 9.52437*((obsfile->obsdata[j].one_obs_data[prn_pos[p]].p2 - obsfile->obsdata[j].one_obs_data[prn_pos[p]].p1 + dcb_sta + dcb_sat)*w + (1.0 - w)*(vtec0 / 9.52437*mf0 +
				 ((obsfile->obsdata[j].one_obs_data[prn_pos[p]].l1*c / f1 - obsfile->obsdata[j].one_obs_data[prn_pos[p]].l2*c / f2) - (obsfile->obsdata[j - 1].one_obs_data[prn_pos[p + 1]].l1*c / f1 - obsfile->obsdata[j - 1].one_obs_data[prn_pos[p + 1]].l2*c / f2)))) / (*mf);
		 }
		 if (vtec1<0.0){
			 qq = 0;
			 continue;
		 }
		 num++;
		 qq++;
		 mf0 = *mf;
		 vtec0 = vtec1;
	 }
	 delete sate_coor;
	 delete pep;
	 delete mf;
	 return vtec1;
 }
 //����VTEC
 void vtec(pobs obsfile, int n, psp3 sp3file, pio ionfile, pdcb dcbfile, psh sh_file, pvt pt, string a0_path){
	 one_result re;
	 BLH ipp;
	 double mf;
	 double vtec_obs,vtec_sh,vtec_ion,a0;
	 XYZ sta_coor;//��վwgs84����
	 BLH sta_coor_blh;//��վ�ʹ��̵�Ĵ������
	 XYZ sate_coor;//��������
	 ENUPOLAR pep;
	 JULIANDAY jul;
	 //double legend[250];
	 double mjd;
	 double pole_lat = 80.33, pole_lon = -72.67;
	 int num = 0;
	 //vector<double>b;
	 //MatrixXd  mb(400, 1), atb(256, 1), ma(400, 256), ata(256, 256), mx(256, 1), mx1(256, 1), resdul(400, 1),inv_ata(256,256);
	 MatrixXd  atb((O + 1)*(O + 1), 1), ata((O + 1)*(O + 1), (O + 1)*(O + 1)), mx((O + 1)*(O + 1), 1), mx1((O + 1)*(O + 1), 1), inv_ata((O + 1)*(O + 1), (O + 1)*(O + 1));
	 //double A[][256] = { 0.0 };
	 //double b[] = { 0.0 };
	 //vector<double *>A;
	 vector<one_legendre>A;
	 vector<double>b;
	 //double a[256];
	 one_legendre a;
	 //double ATA[256][256] = { 0.0 };
	 //��̬�����ڴ�
	/* int m = 400, n = 256;
	 double **A = new double*[m]; //������  
	 for (int i = 0; i < m; i++)
		 A[i] = new double[n]; //������  
	 double **ATA = new double*[n]; //������  
	 for (int i = 0; i < n; i++)
		 ATA[i] = new double[n]; //������  */
	 //int n = 256;
	 double **ATA = new double*[9]; //������  
	 for (int i = 0; i < 9; i++)
		 ATA[i] = new double[9]; //������ 
	 double **ctable = new double*[25]; //������  
	 for (int i = 0; i < 25; i++)
		 ctable[i] = new double[17]; //������ 
	 read_a0(a0_path, ctable);
	 double ATb[9] = { 0.0 };
	 double x[9] = { 0.0 };
	 //cout << "��ʼѭ��" << endl;
	 //UTC day_first;//һ���0ʱ,30sһ��������ӣ���ͬ��վ��n���۲����ݵĹ۲�ʱ�䲻һ��
	 //day_first = obsfile[0].obsdata[0].utime_o;
	 //day_first.hour = 0; day_first.minute = 0; day_first.second = 0.0;
	 ofstream outfile("E:\\GNSS\\����\\vtec.txt", ios::out);
	 ofstream outfile1("E:\\GNSS\\����\\arry.txt", ios::out);
	 ofstream outfile2("E:\\GNSS\\����\\ata.txt", ios::out);
	 ofstream outfile3("E:\\GNSS\\����\\atb.txt", ios::out);
	 ofstream outfile4("E:\\GNSS\\����\\inv_ata.txt", ios::out);
	// outfile.precision(4);
	// outfile1.precision(4);
	 if (!outfile)
	 {
		 cerr << "�ļ�����ʧ��" << endl;
		 abort();
	 }
	 for (int i = 0; i < 30; i++)//һ��֮�еĹ۲������30sһ�Σ���2880��
	 {
		 outfile.setf(ios_base::fixed, ios_base::floatfield);
		 outfile.precision(12);
		 outfile1.setf(ios_base::fixed, ios_base::floatfield);
		 outfile1.precision(12);
		 outfile2.setf(ios_base::fixed, ios_base::floatfield);
		 outfile2.precision(8);
		 outfile3.setf(ios_base::fixed, ios_base::floatfield);
		 outfile3.precision(8);
		 outfile4.setf(ios_base::fixed, ios_base::floatfield);
		 outfile4.precision(8);
		 //outfile5.setf(ios_base::fixed, ios_base::floatfield);
		 //outfile5.precision(8);
		 outfile << setw(4) << obsfile[0].obsdata[i].utime_o.year << setw(3) << obsfile[0].obsdata[i].utime_o.month << setw(3) << obsfile[0].obsdata[i].utime_o.day <<
			 setw(3) << obsfile[0].obsdata[i].utime_o.hour << setw(3) << obsfile[0].obsdata[i].utime_o.minute << setw(15) << obsfile[0].obsdata[i].utime_o.second << endl;
		 re.rtime = obsfile[0].obsdata[i].utime_o;
		 //re.rtime = day_first;
		 //memset(b, 0.0, sizeof(b));
		 //memset(A, 0.0, sizeof(A));
		 //double A[][256] = { 0.0 };
		 //double b[] = { 0.0 };
		 A.clear();
		 b.clear();
		 re.residul.clear();
		// b.clear();
		 num = 0;
		 //cout << "ʱ�䣺" << day_first.hour << ":" << day_first.minute << ":" << day_first.second << endl;
		 cout << "ʱ�䣺" << obsfile[0].obsdata[i].utime_o.hour << ":" << obsfile[0].obsdata[i].utime_o.minute << ":" << obsfile[0].obsdata[i].utime_o.second << endl;
		 for (int j = 0; j < n; j++)//n����վ���������в�վ��ͬһ��ʱ��۲���������һ��������������г����ϵ��
		 {
			 if (fabs(deltjulianday(&obsfile[0].obsdata[i].utime_o, &obsfile[j].obsdata[i].utime_o)) > 0.1)
			 {
				 cout << "                                                                      ʱ���в���" << endl;
				 continue; }//30sһ���۲�����һ����˵o�ļ���30s����һ�ι۲�����
			 sta_coor = obsfile[j].obsheaddata.approx_coordinate;
			 xyztoblh(&sta_coor, &sta_coor_blh);
			 //cout << "��վ����" << obsfile[j].obsheaddata.station << endl;
			 //cout << "ʱ�䣺" << obsfile[j].obsdata[i].utime_o.hour << ":" << obsfile[j].obsdata[i].utime_o.minute << ":" << obsfile[0].obsdata[i].utime_o.second << endl;
			 for (int k = 0; k < obsfile[j].obsdata[i].prn_list.size(); k++)//һ����վ��һ��ʱ��۲⵽����������
			 {
				 if (!cal_sp3_sate_coor(obsfile[j].obsdata[i].prn_list[k], &obsfile[j].obsdata[i].utime_o, sp3file, &sate_coor))//����������
				 {
					 continue;
				 }//�����������ʧ��
				 sate_azi_ele(&sta_coor, &sate_coor, &pep);//�����Ƿ�λ�Ǻ͸߶Ƚ�
				 if (pep.elevation*180.0 / PI < 15.0){ continue; }//���Ǹ߶Ƚ�Ҫ����15��
				 ipp_pos(&sta_coor_blh, &pep, &ipp, &mf);//�󴩴̵�γ�Ⱦ���,��λ����
				 //vtec_obs = obs_vtec(&obsfile[j], obsfile[j].obsdata[i].utime_o, i, obsfile[j].obsdata[i].prn_list[k], sp3file, dcbfile, &ipp);
				 vtec_ion = ionex_vtec(ionfile, obsfile[j].obsdata[i].utime_o, ipp.latitude*180.0 / PI, ipp.longitude *180.0 / PI) / 10.0;
				 vtec_sh = sh_vtec(sh_file, obsfile[j].obsdata[i].utime_o, ipp.latitude, ipp.longitude,15);
				 if (vtec_sh == 0.0){ cout << "����sh_vtecʧ��" << endl; continue; }
				 utctojulianday(&obsfile[j].obsdata[i].utime_o, &jul);//�Ѽ��㴩�̵�ʱ��utcת��������
				 mjd = jul.daynum + (jul.secondnum + jul.secondfrac) / 86400.0 - 2400000.5;//������ת��Լ��������
				 g2m(&ipp, mjd, pole_lat*PI / 180.0, pole_lon*PI / 180.0);//���̵��ɴ������ת���չ̵ش�����
				 legendre(ipp.latitude, ipp.longitude, a.le,O);//�ɴ��̵�������õ¶���ʽ��ֵ
				 //cout << obsfile[j].obsdata[i].prn_list[k] << "vtec:" << vtec_sh << "  " << vtec_ion << endl;
				 a0 = c_a0(ctable, obsfile[j].obsdata[i].utime_o, ipp.latitude, ipp.longitude);//a0
				 //cout << "a0:" << a0 << endl;
				 //if (vtec_ion < 0.0 || vtec_sh < 0.0){ continue; }
				 //for (int l = 0; l < 255; l++){
				//	 outfile1 << std::right << setw(10) << A[num][l]<<",";
					//  }
				 //outfile1 << std::right << setw(10) << A[num][255] << ";" << endl;
				 //b[num] = vtec_sh;
				 //cout << obsfile[j].obsdata[i].prn_list[k] << " vtec_sh-a0:" << vtec_sh - a0 << endl;
				 b.push_back(vtec_sh-a0);
				 A.push_back(a);
				 outfile << std::right << setw(20) << ipp.latitude << std::right << setw(20) << ipp.longitude << std::right << setw(20) << vtec_ion << std::right << setw(20) << vtec_sh <<std::right << 
					 setw(20) << vtec_ion - a0 << std::right << setw(20) << vtec_sh - a0 << std::right << setw(20) << obsfile[j].obsheaddata.station << std::right << setw(20) << obsfile[j].obsdata[i].prn_list[k] << endl;
				 //b[num] = vtec_ion;
				 //mb(num, 0) = vtec_sh;
				 num++;
				 //cout << obsfile[j].obsdata[i].prn_list[k] << "vtec:" << vtec_sh <<"  "<< vtec_ion << endl;
				 //if (num >= 256){ break; }
				 //if (i == 1){
				//	cout << obsfile[j].obsdata[i].prn_list[k] << "  vtec:" << vtec_sh <<"  "<< vtec_ion << endl; 
				 //}
			 }
			// if (num >= 256){ break; }
		 }
		 cout << "�۲������" << num << endl;
		 //outfile1 << endl;
		 //�������
		if (num < (O+1)*(O+1)){ continue; }
		MatrixXd  mb(num, 1), ma(num, (O + 1)*(O + 1)), resdul(num, 1);
		for (int i = 0; i < num; i++){
			mb(i, 0) = b[i];
			for (int j = 0; j < (O + 1)*(O + 1); j++){
				 ma(i,j)= A[i].le[j];//���鸳ֵ����Ӧ�ľ������ھ�������
				 //cout <<"ma"<< ma(i, j) << endl;
				 outfile1 << std::right << setw(20) << A[i].le[j];
				 if (j == (O + 1)*(O + 1)-1){ outfile1 << ";" << endl; }
				 else{ outfile1 << ","; }
			 }
		 }
		outfile1 <<endl;
		 ata = ma.transpose()*ma;
		 atb = ma.transpose()*mb;
		 inv_ata = ata.inverse();
		 mx1 = (ata.inverse())*atb;//ֱ�ӷ���������ⷽ��
		 for (int i = 0; i < (O + 1)*(O + 1); i++)
		 {
			 ATb[i] = atb(i, 0);//����ֵ����Ӧ������
			 outfile3 << std::right << setw(15) << ATb[i];
		 }
		 for (int i = 0; i < (O + 1)*(O + 1); i++)
		 {
			 for (int j = 0; j < (O + 1)*(O + 1); j++)
			 {
				  ATA[i][j] = ata(i, j);//����ֵ����Ӧ������
				  outfile2 << std::right << setw(15) << ATA[i][j];
				  if (j == (O + 1)*(O + 1)-1){ outfile2 << ";" << endl; }
				  else{ outfile2 << ","; }
				  outfile4 << std::right << setw(25) << inv_ata(i,j);
				  if (j == (O + 1)*(O + 1)-1){ outfile4 << ";" << endl; }
				  else{ outfile4 << ","; }
			 }
		 }
		 outfile2 << endl;
		 outfile3 << endl;
		 outfile4 << endl;
		 chol_eq(ATA, ATb, x,(O+1)*(O+1));//cholesky�ֽⷨ�ⷽ��
		 for (int i = 0; i < (O + 1)*(O + 1); i++){
			 mx(i, 0)=x[i];//�������֮�󷴸�ֵ����Ӧ�ľ���
			 re.coe[i] = x[i];//��
			 re.coe1[i] = mx1(i, 0);
		 }
		 resdul = ma*mx - mb;//����������в�
		 for (int a = 0; a < num; a++){
			 re.residul.push_back(resdul(a, 0));//�в�
		 }
		 pt->allresult.push_back(re);
	 }

	 outfile.close();
	 outfile1.close();
	 outfile2.close();
	 outfile3.close();
	 outfile4.close();
	 for (int i = 0; i < 9; i++)
		 delete[] ATA[i];
	 delete[] ATA;
	 for (int i = 0; i < 25; i++)
		 delete[] ctable[i];
	 delete[] ctable;
	 cout << "�������" << endl;
	 cout << "���������" << pt->allresult.size() << endl;
 }

