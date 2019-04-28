
function epfscore = epf(img1,img2)

  
    img1=double(img1);
    img2=double(img2);
    
    h= fspecial('laplacian');
    
    del_img1=imfilter(img1,h);
    del_img2=imfilter(img2,h);
    
    mu_del_img1=mean2(del_img1);%mean of img1
    mu_del_img2=mean2(del_img2);%mean of img 2

    num1=del_img1-mu_del_img1;
    num2=del_img2-mu_del_img2;
    
    n1=num1.*num2;
    
    den1=(num1.*num1);
    d1=sum(den1(:));
    den2=(num2.*num2);
    d2=sum(den2(:));
    
    n=sum(n1(:));
    d=sqrt(d1*d2);
    
    epfscore=n/d;
    
end
