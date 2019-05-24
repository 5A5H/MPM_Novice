#include "sms.h"

double SMSDot(double a[],double b[],int n)
 {double s=0e0;int i;
  for(i=0;i<n;i++)s+=a[i]*b[i];
  return(s);
 }


double SMSSum(double a[],int n)
 {double s=0e0;int i;
  for(i=0;i<n;i++)s+=a[i];
  return(s);
 }


double SMSKDelta(int i,int j)
{if(i==j){return 1e0;}else{return 0e0;};}


double SMSDeltaPart(double *a,int i,int j,int k)
{div_t d=div(i,j);
  if(d.rem || d.quot>k){return 0e0;}else{return a[d.quot-1];};
}

void SMSMove(double a[],double b[],int n)
{int i;
 for(i=0;i<n;i++)b[i]=a[i];
}

void SMSZero(double a[],int n)
{int i;
 for(i=0;i<n;i++)a[i]=0e0;
}
