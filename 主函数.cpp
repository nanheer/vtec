#include"����.h"
#include"stdafx.h"
#include<iostream>
#include<string>
#include"stdlib.h"
#include <io.h>
#include <Eigen/Dense>
using Eigen::MatrixXd;
int _tmain(int argc, _TCHAR* argv[])
{
	/*JULIANDAY jul;
	UTC day_first;//һ���0ʱ,30sһ��������ӣ���ͬ��վ��n���۲����ݵĹ۲�ʱ�䲻һ��
	day_first.year = 2019; day_first.month = 8; day_first.day = 16;
	day_first.hour = 0; day_first.minute = 0; day_first.second = 0.0;
	for (int i = 0; i < 10; i++){
		utctojulianday(&day_first, &jul);
		jul.secondnum += 30;
		if (jul.secondnum + jul.secondfrac >= 86400.0){ jul.daynum += 1; jul.secondnum -= 86400; }
		juliandaytoutc(&jul, &day_first);\
		cout << "ʱ�䣺" << day_first.hour << ":" << day_first.minute << ":" << day_first.second << endl;
	}*/
	/*MatrixXd  mb(10, 1), atb(4, 1), ma(10, 4), ata(4, 4), mx(4, 1), mx1(4, 1), resdul(10, 1);
	double A[10][4] = { 2, -1, 3, 1, 4, 2, 5, 1, 2, 0, 2, 1, 1, 1, 1, 1, 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	double b[10] = { 2, 5, 7, 3, -7, 0, 0, 0, 0, 0 };
	double ATA[4][4] = { 0.0 };
	double ATb[4] = { 0.0 };
	double x[4] = { 0.0 };
	//MatrixXd  mb(10, 1), ma(10, 4),x(4,1);
	for (int i = 0; i < 10; i++){
		mb(i, 0) = b[i];
		for (int j = 0; j < 4; j++){
			ma(i, j) = A[i][j];//���鸳ֵ����Ӧ�ľ������ھ�������
			//cout <<"ma"<< ma(i, j) << endl;
		}
	}
	ata = ma.transpose()*ma;
	atb = ma.transpose()*mb;
	mx1 = (ata.inverse())*atb;//ֱ�ӷ���������ⷽ��
	for (int i = 0; i < 4; i++)
	{
		ATb[i] = atb(i, 0);//����ֵ����Ӧ������
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			ATA[i][j] = ata(i, j);//����ֵ����Ӧ������
		}
	}
	chol_eq(ATA, ATb, x);//cholesky�ֽⷨ�ⷽ��
	for (int i = 0; i < 4; i++){
		cout << "cholesky��⣺" << x[i] << endl;
		cout << "����ֱ����⣺" << mx1(i,0) << endl;
		mx(i, 0) = x[i];//�������֮�󷴸�ֵ����Ӧ�ľ���;
	}
	resdul = ma*mx - mb;//����������в�
	for (int i = 0; i < 10; i++){
		cout << "�в" << resdul(i, 0)<< endl;
	}*/
	pobs po = new obs;
	obs obs_data;
	psp3 pn = new all_sate_ephem;
	pio ion = new ionex;
	pdcb pd = new dcb;
	psh ps = new sh_file;
	pvt pt = new result;
	vector<obs> sta_obs;
	//string str_file="E:/GNSS/�����/CODE/2001";
	//string all_str[150];
	obs obs_value[40];
	//obs obs_value;
	//dir(str_file,obs_value);
	//put_file(str_file, all_str, 0);
	//������ȡ����o�ļ�
	string path[40] = 
	{ 
		"E:/GNSS/����/o�ļ�/bako0240.17o", "E:/GNSS/����/o�ļ�/bamf0240.17o", "E:/GNSS/����/o�ļ�/barh0240.17o", "E:/GNSS/����/o�ļ�/bcov0240.17o",
		"E:/GNSS/����/o�ļ�/cebr0240.17o",
		"E:/GNSS/����/o�ļ�/chu20240.17o", "E:/GNSS/����/o�ļ�/chum0240.17o", "E:/GNSS/����/o�ļ�/chur0240.17o", "E:/GNSS/����/o�ļ�/chwk0240.17o",
		"E:/GNSS/����/o�ļ�/cibg0240.17o", "E:/GNSS/����/o�ļ�/cit10240.17o", "E:/GNSS/����/o�ļ�/ckis0240.17o", "E:/GNSS/����/o�ļ�/clrs0240.17o",
		"E:/GNSS/����/o�ļ�/brst0240.17o",
		"E:/GNSS/����/o�ļ�/brux0240.17o", "E:/GNSS/����/o�ļ�/bucu0240.17o",
		"E:/GNSS/����/o�ļ�/budp0240.17o", "E:/GNSS/����/o�ļ�/bzrg0240.17o", "E:/GNSS/����/o�ļ�/cags0240.17o", 
		"E:/GNSS/����/o�ļ�/cata0240.17o", "E:/GNSS/����/o�ļ�/ccj20240.17o", "E:/GNSS/����/o�ļ�/chan0240.17o",
		"E:/GNSS/����/o�ļ�/chil0240.17o", "E:/GNSS/����/o�ļ�/chpg0240.17o", "E:/GNSS/����/o�ļ�/chpi0240.17o", "E:/GNSS/����/o�ļ�/chti0240.17o",

		"E:/GNSS/����/o�ļ�/albh0240.17o", "E:/GNSS/����/o�ļ�/alg20240.17o", "E:/GNSS/����/o�ļ�/arev0240.17o", "E:/GNSS/����/o�ļ�/bake0240.17o",
		"E:/GNSS/����/o�ļ�/bik00240.17o", "E:/GNSS/����/o�ļ�/cas10240.17o", "E:/GNSS/����/o�ļ�/cebr0240.17o", "E:/GNSS/����/o�ļ�/dlf10240.17o",
		"E:/GNSS/����/o�ļ�/drao0240.17o", "E:/GNSS/����/o�ļ�/ganp0240.17o", "E:/GNSS/����/o�ļ�/shao0240.17o",
		"E:/GNSS/����/o�ļ�/auck0240.17o", "E:/GNSS/����/o�ļ�/badg0240.17o", "E:/GNSS/����/o�ļ�/baie0240.17o"
		/*
		"E:/GNSS/����/o�ļ�/bhr30240.17o", "E:/GNSS/����/o�ļ�/bhr40240.17o", "E:/GNSS/����/o�ļ�/bilb0240.17o", "E:/GNSS/����/o�ļ�/bill0240.17o",
		"E:/GNSS/����/o�ļ�/bin10240.17o", "E:/GNSS/����/o�ļ�/bjnm0240.17o", "E:/GNSS/����/o�ļ�/blyt0240.17o", "E:/GNSS/����/o�ļ�/bogi0240.17o",
		"E:/GNSS/����/o�ļ�/bogt0240.17o", "E:/GNSS/����/o�ļ�/bor10240.17o", "E:/GNSS/����/o�ļ�/braz0240.17o", "E:/GNSS/����/o�ļ�/brew0240.17o",
		"E:/GNSS/����/o�ļ�/brft0240.17o"*/                                
	};
	
	for (int i = 0; i < 40; i++)
	{
		readofile_vtec(path[i], &obs_value[i]);
	}
	string strn = "E:\\GNSS\\����\\igs19332.sp3";
	string stro = "E:\\GNSS\\����\\cebr0240.17o";//cas1
	string stri = "E:\\GNSS\\����\\codg0240.17i";
	string strd = "E:\\GNSS\\����\\CAS0MGXRAP_20170240000_01D_01D_DCB.BSX";
	string strs = "E:\\GNSS\\�����\\CODE\\2017\\COD19332.ION";
	string stra = "E:\\GNSS\\����\\a0.dat";
	//readofile_vtec(stro, po);
	//sta_coor =po->obsheaddata.approx_coordinate;
	//xyztoblh(&sta_coor, &sta_coor_blh);
	//cout << "��վ�������꣺"<<po->obsheaddata.station<<"��" << sta_coor_blh.latitude*180.0 / PI << "," << sta_coor_blh.longitude*180.0 / PI << ")" <<  endl;
	//return 0;
	//readofile_vtec(stro, po);
	readsp3file(strn, pn);
	read_ionex(stri, ion);
	read_dcb(strd, pd);
	read_sh(strs, ps);
	vtec(obs_value,40, pn,ion,pd,ps,pt,stra);
	putresult(pt);
	delete po;
	//delete obs_value;
	delete pn;
	delete ion;
	delete pd;
	delete ps;
	delete pt;
	return 0;
}