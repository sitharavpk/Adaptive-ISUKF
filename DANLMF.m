
function output = DANLMF(damrf,sigma)

sigma=sigma(sigma>0);
sigma_n=50;

input1=damrf;
[m n]=size(input1);

t=9;
f=7;

% Memory for the output
 output=zeros(m,n);

% Replicate the boundaries of the input image
 input2 = padarray(input1,[f f],'symmetric');

% Used kernel
 kernel = make_kernel(f);
 kernel = kernel / sum(sum(kernel));

 for i=1:m
 for j=1:n
                 
         i1 = i+ f;
         j1 = j+ f;
                
         W1= input2(i1-f:i1+f , j1-f:j1+f);
         
         wmax=0; 
         average=0;
         sweight=0;
         
         rmin = max(i1-t,f+1);
         rmax = min(i1+t,m+f);
         smin = max(j1-t,f+1);
         smax = min(j1+t,n+f);
         
         for r=rmin:1:rmax
         for s=smin:1:smax
                                               
                if(r==i1 && s==j1) continue; end;
                                
                W2= input2(r-f:r+f , s-f:s+f);                
                
                d = sum(sum(kernel.*(W1-W2).*(W1-W2)));
                
                w=exp(-d/50);
                
                if w>wmax                
                    wmax=w;                   
                end
                
                sweight = sweight + w;
                average = average + w*input2(r,s);                                  
         end 
         end
             
        average = average + wmax*input2(i1,j1);
        sweight = sweight + wmax;
                   
        if sweight > 0
            output(i,j) = average / sweight;
        else
            output(i,j) = input1(i,j);
        end                
 end
 end
