#ifndef FUNCTION_H
#define FUNCTION_H
#include<string>
#include"�ṹ��.h"
#include"����.h"
//#include"matrix.h"
using namespace std;
//��sp3�ļ�
void readsp3file(string strn, psp3 sp3file);
//����vtec����ķ�ʽ��ȡo�ļ�
void readofile_vtec(string stro, pobs obsfile);
//��ȡionex�ļ�
void read_ionex(string stri, pio ionfile);
//��ȡdcb�ļ�
void read_dcb(string strd, pdcb dcbfile);
//GPST��UTC
void gpsttoutc(pgpst gt, ptc ut);
//JULIANDAY��BDT
void juliandaytobdt(pjulian ju, pbdt bt);
//UTC��BDT
void utctobdt(ptc ut, pbdt bt);
//UTC��GPST
void utctogpst(ptc ut,pgpst gt);

//UTC��JULIANDAY
void utctojulianday(ptc ut,pjulian ju);

//JULIANDAY��UTC
void juliandaytoutc(pjulian ju,ptc ut);

//GPST��JULIANDAY
void gpsttojulianday(pgpst gt, pjulian ju);

//JULIANDAY��GPST
void juliandaytogpst(pjulian ju, pgpst gt);
//�������ղ�ֵ
double deltjulianday(ptc u1, ptc u2);
//�����ձ任
void transjulian(pjulian ju1, double * dt, pjulian ju2);
//

//XYZ��BLH
void xyztoblh(pxyz px,pblh pb);
//BLH��XYZ
void blhtoxyz(pblh pb, pxyz px);
//XYZ��ENU
void xyztoenu(pxyz pxcenter, pxyz px, penu pe);
//ENU��XZY
void enutoxyz(pxyz pxcenter,penu pe, pxyz px);
//ENU��ENUPOLAR
void enutoenupolar(penu pe, penupolar pep);
//�����Ƿ�λ�Ǻ͸߶Ƚ�
void sate_azi_ele(pxyz pxcenter, pxyz px, penupolar pep);
//�������ӳٸ���
void tropo(double mjd,pxyz px1, pxyz px2,double * dx);
//����������
void get_nominal_metedata(double mjd, double lat, double lon,double dhgt, double* pres, double* temp, double* rhumi, double* undu);
//���Сֵ
int min(int a, int b);
//���ջ����߸߸���
void antena_height_correction(pxyz px1, pxyz px2, penu pe, double*p1, double*p2, double*pp1, double*pp2);
//��������������������
bool cal_sate_coor(string prn, ptc ut, psp3 sp3file, pxyz coor);
//��ƽ����Ǽ���ƫ�����
double calE(double* M, double* e);
//��ƫ�������������
double calf(double* E, double* e);
//��ת����
//void rotamatrix(string str, double seta,Matrix r);
//������
void putresult(pvt pt);
//vtec����
void vtec(pobs obsfile, int n, psp3 sp3file, pio ionfile, pdcb dcbfile, psh sh_file, pvt pt, string a0_path);
//��ionex�ļ���ֵ����������vtec
double ionex_vtec(pio ionfile, UTC vtime, double lat, double lon);//lat��lon�Ǵ�ؾ�γ��
//��dcb�ļ����ҵ���Ӧ��dcbֵ
bool find_dcb(bool b, string prn, string station, pdcb dcbfile, double* dcb_val);
//bool find_dcb(bool b, string prn, string station, pio ionfile, double* dcb_val);
//��ȡ�������г�����ļ�
void read_sh(string path, psh sh_file);
//�ɵ������г�����ļ���������vtec
double sh_vtec(psh sh_file, UTC vtime, double lat, double lon,int n);
//����黯���õº�����ֵ
void legendre(double lat, double**pg, double *pz,int n);
//����黯���õ¶���ʽ��г��
void geo_legendre(double lat, double**pg, int n);
//����黯���õ¶���ʽ��г��
void zone_legendre(double lat, double *pz,int n);
//������굽�չ̵ش�����
void g2m(pblh pb1, double pole_lat, double pole_lon);
void Transposition(double**L, int n);//����ת��
void Multi(double**A, double**B, int n);//�������
void chol_eq(double**A, double* b, double* x, int n);//Cholesky�ֽⷽ���ⷽ����
void chol_rf(double**A, double**L, double**d, int n);//Cholesky�ֽⷽ��
//�۲�ֵ����vtec
double obs_vtec(pobs obsfile, UTC obs_t, int n, string prn, psp3 sp3file, pdcb dcbfile, pblh ipp);
void put_file(string path, string* all_file, int layer);
void dir(string path, string* all_file);
void read_a0(string path, double** ctable);//��ȡ����ϵ����
double c_a0(double**ctable, UTC ut, double lat, double lon);
#endif