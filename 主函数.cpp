#include"函数.h"
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
	UTC day_first;//一天的0时,30s一个间隔增加，不同测站第n个观测数据的观测时间不一致
	day_first.year = 2019; day_first.month = 8; day_first.day = 16;
	day_first.hour = 0; day_first.minute = 0; day_first.second = 0.0;
	for (int i = 0; i < 10; i++){
		utctojulianday(&day_first, &jul);
		jul.secondnum += 30;
		if (jul.secondnum + jul.secondfrac >= 86400.0){ jul.daynum += 1; jul.secondnum -= 86400; }
		juliandaytoutc(&jul, &day_first);\
		cout << "时间：" << day_first.hour << ":" << day_first.minute << ":" << day_first.second << endl;
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
			ma(i, j) = A[i][j];//数组赋值给对应的矩阵用于矩阵运算
			//cout <<"ma"<< ma(i, j) << endl;
		}
	}
	ata = ma.transpose()*ma;
	atb = ma.transpose()*mb;
	mx1 = (ata.inverse())*atb;//直接方程组求逆解方程
	for (int i = 0; i < 4; i++)
	{
		ATb[i] = atb(i, 0);//矩阵赋值给对应的数组
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			ATA[i][j] = ata(i, j);//矩阵赋值给对应的数组
		}
	}
	chol_eq(ATA, ATb, x);//cholesky分解法解方程
	for (int i = 0; i < 4; i++){
		cout << "cholesky求解：" << x[i] << endl;
		cout << "矩阵直接求解：" << mx1(i,0) << endl;
		mx(i, 0) = x[i];//数组求解之后反赋值给对应的矩阵;
	}
	resdul = ma*mx - mb;//矩阵运算求残差
	for (int i = 0; i < 10; i++){
		cout << "残差：" << resdul(i, 0)<< endl;
	}*/
	pobs po = new obs;
	obs obs_data;
	psp3 pn = new all_sate_ephem;
	pio ion = new ionex;
	pdcb pd = new dcb;
	psh ps = new sh_file;
	pvt pt = new result;
	vector<obs> sta_obs;
	//string str_file="E:/GNSS/电离层/CODE/2001";
	//string all_str[150];
	obs obs_value[40];
	//obs obs_value;
	//dir(str_file,obs_value);
	//put_file(str_file, all_str, 0);
	//遍历读取所有o文件
	string path[40] = 
	{ 
		"E:/GNSS/数据/o文件/bako0240.17o", "E:/GNSS/数据/o文件/bamf0240.17o", "E:/GNSS/数据/o文件/barh0240.17o", "E:/GNSS/数据/o文件/bcov0240.17o",
		"E:/GNSS/数据/o文件/cebr0240.17o",
		"E:/GNSS/数据/o文件/chu20240.17o", "E:/GNSS/数据/o文件/chum0240.17o", "E:/GNSS/数据/o文件/chur0240.17o", "E:/GNSS/数据/o文件/chwk0240.17o",
		"E:/GNSS/数据/o文件/cibg0240.17o", "E:/GNSS/数据/o文件/cit10240.17o", "E:/GNSS/数据/o文件/ckis0240.17o", "E:/GNSS/数据/o文件/clrs0240.17o",
		"E:/GNSS/数据/o文件/brst0240.17o",
		"E:/GNSS/数据/o文件/brux0240.17o", "E:/GNSS/数据/o文件/bucu0240.17o",
		"E:/GNSS/数据/o文件/budp0240.17o", "E:/GNSS/数据/o文件/bzrg0240.17o", "E:/GNSS/数据/o文件/cags0240.17o", 
		"E:/GNSS/数据/o文件/cata0240.17o", "E:/GNSS/数据/o文件/ccj20240.17o", "E:/GNSS/数据/o文件/chan0240.17o",
		"E:/GNSS/数据/o文件/chil0240.17o", "E:/GNSS/数据/o文件/chpg0240.17o", "E:/GNSS/数据/o文件/chpi0240.17o", "E:/GNSS/数据/o文件/chti0240.17o",

		"E:/GNSS/数据/o文件/albh0240.17o", "E:/GNSS/数据/o文件/alg20240.17o", "E:/GNSS/数据/o文件/arev0240.17o", "E:/GNSS/数据/o文件/bake0240.17o",
		"E:/GNSS/数据/o文件/bik00240.17o", "E:/GNSS/数据/o文件/cas10240.17o", "E:/GNSS/数据/o文件/cebr0240.17o", "E:/GNSS/数据/o文件/dlf10240.17o",
		"E:/GNSS/数据/o文件/drao0240.17o", "E:/GNSS/数据/o文件/ganp0240.17o", "E:/GNSS/数据/o文件/shao0240.17o",
		"E:/GNSS/数据/o文件/auck0240.17o", "E:/GNSS/数据/o文件/badg0240.17o", "E:/GNSS/数据/o文件/baie0240.17o"
		/*
		"E:/GNSS/数据/o文件/bhr30240.17o", "E:/GNSS/数据/o文件/bhr40240.17o", "E:/GNSS/数据/o文件/bilb0240.17o", "E:/GNSS/数据/o文件/bill0240.17o",
		"E:/GNSS/数据/o文件/bin10240.17o", "E:/GNSS/数据/o文件/bjnm0240.17o", "E:/GNSS/数据/o文件/blyt0240.17o", "E:/GNSS/数据/o文件/bogi0240.17o",
		"E:/GNSS/数据/o文件/bogt0240.17o", "E:/GNSS/数据/o文件/bor10240.17o", "E:/GNSS/数据/o文件/braz0240.17o", "E:/GNSS/数据/o文件/brew0240.17o",
		"E:/GNSS/数据/o文件/brft0240.17o"*/                                
	};
	
	for (int i = 0; i < 40; i++)
	{
		readofile_vtec(path[i], &obs_value[i]);
	}
	string strn = "E:\\GNSS\\数据\\igs19332.sp3";
	string stro = "E:\\GNSS\\数据\\cebr0240.17o";//cas1
	string stri = "E:\\GNSS\\数据\\codg0240.17i";
	string strd = "E:\\GNSS\\数据\\CAS0MGXRAP_20170240000_01D_01D_DCB.BSX";
	string strs = "E:\\GNSS\\电离层\\CODE\\2017\\COD19332.ION";
	string stra = "E:\\GNSS\\数据\\a0.dat";
	//readofile_vtec(stro, po);
	//sta_coor =po->obsheaddata.approx_coordinate;
	//xyztoblh(&sta_coor, &sta_coor_blh);
	//cout << "测站地理坐标："<<po->obsheaddata.station<<"（" << sta_coor_blh.latitude*180.0 / PI << "," << sta_coor_blh.longitude*180.0 / PI << ")" <<  endl;
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