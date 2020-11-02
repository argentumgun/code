/* This code accompanies
 *   The Lattice Boltzmann Method: Principles and Practice
 *   T. Krüger, H. Kusumaatmaja, A. Kuzmin, O. Shardt, G. Silva, E.M. Viggen
 *   ISBN 978-3-319-44649-3 (Electronic) 
 *        978-3-319-44647-9 (Print)
 *   http://www.springer.com/978-3-319-44647-9
 *
 * This code is provided under the MIT license. See LICENSE.txt.
 *
 * Author: Alexandr Kuzmin
 * 
 */
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>

//Domain size
int NY;
int NX;
int NUM;
int radius;

//Other constants，9个速度方向
const int NPOP=9;

//Time steps
int N=64000;
int NOUTPUT=200;

//Fields and populations  流体性质（速度密度）与分布 的 数组
double *f;    //碰撞前分布
double *f2;     //碰撞后分布
double *rho;    //密度
double *ux;     //x方向速度
double *uy;     //y方向速度
int * geometry; //？？？？？？？？？？？？？

//Boundary conditions
double conc_wall=1.0;
std::vector<int> bb_nodes;//边界点
std::vector<char>* dirs;//
int *bottom;
int *top;

//BGK relaxation parameter
double omega=1.0/(0.5+0.5/32.0);
double omega_plus=2.0-omega;
double omega_minus=omega;

//Diffusion parameter
double ce=1.0/3.0;

//Underlying lattice parameters
double weights[]={4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};     //速度分量的权重（从中心点0开始）
double weights_trt[]={0.0,1.0/3.0,1.0/3.0,1.0/3.0,1.0/3.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0};     //？？？
int cx[]={0,1,0,-1,0,1,-1,-1,1};                                                                    //cx cy速度分量（九个速度矢量（x，y））
int cy[]={0,0,1,0,-1,1,1,-1,-1};
int complement[]={0,3,4,1,2,7,8,5,6};                                                               //速度分量的下标
int pxx[]={0, 1, -1, 1, -1, 0, 0, 0, 0};                                                            //
int pxy[]={0, 0, 0, 0, 0, 1, -1, 1, -1};                                                            //

void writedensity(std::string const & fname)
{
    std::string filename = "cylinder/" + fname+".dat";
    std::ofstream fout(filename.c_str());
    fout.precision(10);

    for (int iX=0; iX<NX; ++iX)
    {
        for (int iY=0; iY<NY; iY++)
        {
            int counter=iY*NX+iX;
            fout<<rho[counter]<<" ";//输出格点的密度
        }
        fout<<"\n";
    }
}

void writegeometry(std::string const & fname)
{
    std::string filename = "cylinder/" + fname+".dat";
    std::ofstream fout(filename.c_str());
    fout.precision(10);

    for (int iX=0; iX<NX; ++iX)
    {
        for (int iY=0; iY<NY; iY++)
        {
            int counter=iY*NX+iX;
            fout<<geometry[counter]<<" ";
        }
        fout<<"\n";
    }
}



void init()     //初始化
{
    //Creating arrays               节点个数 * 速度分量个数
    f=new double[NUM*NPOP];         //碰撞前的分布
    f2=new double[NUM*NPOP];        //碰撞后的分布
    
    //Bulk nodes initialization
    double feq;                     //单个节点的
    double geq;
    double sum;
    
    for(int iY=0;iY<NY;iY++)
        for(int iX=0;iX<NX;iX++)    //遍历每个节点
        {
            int  counter=iY*NX+iX;      //计数
            double dense_temp=rho[counter];     //碰撞前密度
            double ux_temp=ux[counter];         //碰撞前速度x方向
            double uy_temp=uy[counter];         //碰撞前速度y方向
            
            for (int k=0; k<NPOP; k++)//遍历每个速度分量
            {
                feq=weights[k]*(dense_temp+3.0*dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp) //计算平衡分布
                                +4.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp
                                                +(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp
                                                +2.0*ux_temp*uy_temp*cx[k]*cy[k]));
                f[counter*NPOP+k]=feq;//将计算出的平衡状态赋值给每个节点的速度分量
            }
            
        }
}



void collide()//碰撞
{
    for (int counter=0;counter<NUM;counter++)//counter 节点编号
        if (geometry[counter]==1)
        {
            rho[counter]=0.0;   //密度初始化为0
            
            int offset=counter*NPOP;    //计算初始下标
     
            double sum=0;
            for(int k=0;k<NPOP;k++)
                sum+=f[offset+k];//节点粒子分布总和
    
            rho[counter]=sum;//碰撞前密度

            double dense_temp=rho[counter];//碰撞前密度
            double ux_temp=ux[counter];//碰撞前x方向速度
            double uy_temp=uy[counter];//碰撞前y方向速度

            //BGK equilibrium
            double feq[NPOP];//声明在BGK算子平衡状态的分布数组
    
            //TRT equilibrium
            double feq_plus[NPOP],feq_minus[NPOP];
            double f_plus[NPOP],f_minus[NPOP];

             //Speeding up the process
            f_plus[0]=f[offset];
            f_plus[1]=0.5*(f[offset+1]+f[offset+3]);
            f_plus[2]=0.5*(f[offset+2]+f[offset+4]);
            f_plus[3]=f_plus[1];
            f_plus[4]=f_plus[2];
            f_plus[5]=0.5*(f[offset+5]+f[offset+7]);
            f_plus[6]=0.5*(f[offset+6]+f[offset+8]);
            f_plus[7]=f_plus[5];
            f_plus[8]=f_plus[6];

            f_minus[0]=0.0;
            f_minus[1]=0.5*(f[offset+1]-f[offset+3]);
            f_minus[2]=0.5*(f[offset+2]-f[offset+4]);
            f_minus[3]=-f_minus[1];
            f_minus[4]=-f_minus[2];
            f_minus[5]=0.5*(f[offset+5]-f[offset+7]);
            f_minus[6]=0.5*(f[offset+6]-f[offset+8]);
            f_minus[7]=-f_minus[5];
            f_minus[8]=-f_minus[6];
          
            feq_minus[0]=0.0;
            feq_minus[1]=weights_trt[1]*dense_temp*(cx[1]*ux_temp+cy[1]*uy_temp);
            feq_minus[2]=weights_trt[2]*dense_temp*(cx[2]*ux_temp+cy[2]*uy_temp);
            feq_minus[3]=-feq_minus[1];
            feq_minus[4]=-feq_minus[2];
            feq_minus[5]=weights_trt[5]*dense_temp*(cx[5]*ux_temp+cy[5]*uy_temp);
            feq_minus[6]=weights_trt[6]*dense_temp*(cx[6]*ux_temp+cy[6]*uy_temp);
            feq_minus[7]=-feq_minus[5];
            feq_minus[8]=-feq_minus[6];
            
            double u_sq=ux_temp*ux_temp+uy_temp*uy_temp;//计算速度的平方
        
            feq_plus[1]=weights_trt[1]*dense_temp*(ce+0.5*(3.0*(cx[1]*ux_temp+cy[1]*uy_temp)*(cx[1]*ux_temp+cy[1]*uy_temp)-u_sq));
            feq_plus[2]=weights_trt[2]*dense_temp*(ce+0.5*(3.0*(cx[2]*ux_temp+cy[2]*uy_temp)*(cx[2]*ux_temp+cy[2]*uy_temp)-u_sq));
            feq_plus[3]=feq_plus[1];
            feq_plus[4]=feq_plus[2];
            feq_plus[5]=weights_trt[5]*dense_temp*(ce+0.5*(3.0*(cx[5]*ux_temp+cy[5]*uy_temp)*(cx[5]*ux_temp+cy[5]*uy_temp)-u_sq));
            feq_plus[6]=weights_trt[6]*dense_temp*(ce+0.5*(3.0*(cx[6]*ux_temp+cy[6]*uy_temp)*(cx[6]*ux_temp+cy[6]*uy_temp)-u_sq));
            feq_plus[7]=feq_plus[5];
            feq_plus[8]=feq_plus[6];
            feq_plus[0]=dense_temp-2.0*(feq_plus[1]+feq_plus[2]+feq_plus[5]+feq_plus[6]);

     
            //Collision operator    通过碰撞算子计算碰撞后分布
            for(int k=0; k < NPOP; k++)
                f2[offset+k]=f[offset+k]-omega_plus*(f_plus[k]-feq_plus[k])-omega_minus*(f_minus[k]-feq_minus[k]);
    }
            
}


void update_bounce_back()   //边界条件
{
    for(int counter=0;counter<bb_nodes.size();counter++)
    {
        for(int k=0;k<dirs[counter].size();k++)
        {
            int dir=dirs[counter][k];
            int counter2=bb_nodes[counter]+cy[dir]*NX+cx[dir];
            f2[bb_nodes[counter]*NPOP+dir]=-f2[counter2*NPOP+complement[dir]]+2*weights[dir]*conc_wall;
        }
    }


}

void initialize_geometry()  //初始化几何边界
{
    NY=129;
    NX=129;
    NUM=NX*NY;
    geometry=new int[NUM];
    rho=new double[NUM];
    ux=new double[NUM];
    uy=new double[NUM];
    radius=40;
    
    for(int iX=0;iX<NX;iX++)
        for(int iY=0;iY<NY;iY++)//遍历每个节点  设置标记
        {
            int counter=iX*NY+iY;
            if ((iX-(NX-1)/2)*(iX-(NX-1)/2)+(iY-(NY-1)/2)*(iY-(NY-1)/2)>radius*radius)   //设置圆柱体范围。内标记为1。外为-1
                geometry[counter]=-1;
            else 
                geometry[counter]=1;
        }
    writegeometry("geometry_before");  //输出到记录初始状态     cylinder/geometry_before.dat  

    for(int iX=0;iX<NX;iX++)
        for(int iY=0;iY<NY;iY++)//遍历每个节点 若与外部相邻则设为边界  
        {
            int counter=iX*NY+iY;
            if (geometry[counter]==1)//暂定此节点在边界内部
            {    
                bool flag=false;// 设置Boolean值 
                for(int k=0;k<NPOP;k++)
                {
                    int counter2=(iX+cx[k])*NY+iY+cy[k];// 遍历此节点的每个相邻节点
                    if (geometry[counter2]==-1)  //相邻节点为外界向左向下或左下方
                        flag=true;
                }
                if (flag)
                     geometry[counter]=0;//设为0
            }
            
        }

    for(int iX=0;iX<NX;iX++)
        for(int iY=0;iY<NY;iY++)
        {
            int counter=iX*NY+iY;
            if (geometry[counter]==0)//若节点为边界
            {    
                bool flag=false;
                for(int k=0;k<NPOP;k++)
                {
                    int counter2=(iX+cx[k])*NY+iY+cy[k];
                    if (geometry[counter2]==1)
                        flag=true;
                }
                if (!flag)
                     geometry[counter]=-1;
            }
            
        }
        
        
    writegeometry("geometry_after");

    for(int counter=0;counter<NUM;counter++)
    {
        if (geometry[counter]==0)  //墙体密度为1 速度为0
        {
            rho[counter]=conc_wall;
            ux[counter]=0.0;
            uy[counter]=0.0;
            bb_nodes.push_back(counter);  //记录进中vector
        }
        else if(geometry[counter]==-1)  //外部 密度设为-1 速度为0
        {
            rho[counter]=-1.0;
            ux[counter]=0.0;
            uy[counter]=0.0;
        }
        else                    //墙体内部 密度速度均为0
        {
            rho[counter]=0.0;
            ux[counter]=0.0;
            uy[counter]=0.0;
        }
    }
    
    
    //Finding directions for BB nodes
    dirs=new std::vector<char>[bb_nodes.size()];
    for(int counter=0;counter<bb_nodes.size();counter++)
    {
        for(int k=1;k<NPOP;k++)
        {
            int counter2=bb_nodes[counter]+cy[k]*NX+cx[k];
            if (geometry[counter2]==1)
                dirs[counter].push_back(k);
        }
    }
}


void finish_simulation()
{
    delete[] geometry;
    delete[] rho;
    delete[] ux;
    delete[] uy;
    delete[] f;
    delete[] f2;
    delete[] dirs;
}


void stream()
{
    for(int counter=0;counter<=NUM;counter++)
        if (geometry[counter]==1)
            for(int iPop=0;iPop<NPOP;iPop++)
            {
                int counter2=counter-cy[iPop]*NX-cx[iPop];
                f[counter*NPOP+iPop]=f2[counter2*NPOP+iPop];
            }

}


int main(int argc, char* argv[])
{
    initialize_geometry();//初始边界
    
    init();//将分布初始为平衡状态
    
    for(int counter=0;counter<=N;counter++)//执行N次循环
    {

        collide();//碰撞
        update_bounce_back();//反弹
        stream();//移动
        
        //Writing files
        if (counter%NOUTPUT==0)//每200次输出一次
        {
            std::cout<<"Counter="<<counter<<"\n";
            std::stringstream filewritedensity;
             
            std::stringstream counterconvert;
            counterconvert<<counter;
            filewritedensity<<std::fixed;

            filewritedensity<<"cylinder"<<std::string(7-counterconvert.str().size(),'0')<<counter;
            
            writedensity(filewritedensity.str());
        }

    }

    finish_simulation();//释放内存
    
    return 0;
}
