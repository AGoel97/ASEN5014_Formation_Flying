function [K, F] = fsf_gains(A, B, C)

%K = place(A,B,...
    %[-0.0010652 -0.0016794 -0.0013118 -.0018838 -.00088127 -.00090389]);
    %[-0.0006734 -0.0013248 -0.0015337 -0.0018912 -0.0019766 -0.0013533]);
    %[-0.0014447 -0.0019018 -0.00079164 -0.0013151 -0.0010307 -0.00079626]);
%K = place(A,B,[-.00079626 -.00079164 -.002*rand(1,4)]);
K = place(A,B,...
    [-0.00110777413157877 -0.00106524700494181 -0.000791640000001095 -0.000903891418517599 -0.000796260000000242 -0.00167939484143928]);
F = inv(C/(-A+B*K)*B);