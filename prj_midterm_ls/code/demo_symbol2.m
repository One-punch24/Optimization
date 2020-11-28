clear all
h = @(c1, c2, c3, c4, c5, t) c5.*c1.^(c2.^((1-c3.^t).^c4));
syms C1 C2 C3 C4 C5 data(i) y x1 x2 x3 x4 x5
C1 = 1/(1+exp(x1));
C2 = 1/(1+exp(x2));
C3 = 1/(1+exp(x3));
C4 = 1+exp(x4)
C5 = x5;
f1 = h(C1, C2, C3, C4, C5, data(i))-y;

df1 = diff(f1, x1)
df2 = diff(f1, x2)
df3 = diff(f1, x3)
df4 = diff(f1, x4)
df5 = diff(f1, x5)

display('-----')
df11 = diff(df1, x1);
df12 = diff(df1, x2);
df13 = diff(df1, x3);
df14 = diff(df1, x4);
df15 = diff(df1, x5);
df21 = diff(df2, x1);
df22 = diff(df2, x2);
df23 = diff(df2, x3);
df24 = diff(df2, x4);
df25 = diff(df2, x5);
df31 = diff(df3, x1);
df32 = diff(df3, x2);
df33 = diff(df3, x3);
df34 = diff(df3, x4);
df35 = diff(df3, x5);
df41 = diff(df4, x1);
df42 = diff(df4, x2);
df43 = diff(df4, x3);
df44 = diff(df4, x4);
df45 = diff(df4, x5);
df51 = diff(df5, x1);
df52 = diff(df5, x2);
df53 = diff(df5, x3);
df54 = diff(df5, x4);
df55 = diff(df5, x5);
display('tmp(1,1) = '+string(df11)+';');
display('tmp(1,2) = '+string(df12)+';');
display('tmp(1,3) = '+string(df13)+';');
display('tmp(1,4) = '+string(df14)+';');
display('tmp(1,5) = '+string(df15)+';');
display('tmp(2,1) = '+string(df21)+';');
display('tmp(2,2) = '+string(df22)+';');
display('tmp(2,3) = '+string(df23)+';');
display('tmp(2,4) = '+string(df24)+';');
display('tmp(2,5) = '+string(df25)+';');
display('tmp(3,1) = '+string(df31)+';');
display('tmp(3,2) = '+string(df32)+';');
display('tmp(3,3) = '+string(df33)+';');
display('tmp(3,4) = '+string(df34)+';');
display('tmp(3,5) = '+string(df35)+';');
display('tmp(4,1) = '+string(df41)+';');
display('tmp(4,2) = '+string(df42)+';');
display('tmp(4,3) = '+string(df43)+';');
display('tmp(4,4) = '+string(df44)+';');
display('tmp(4,5) = '+string(df45)+';');

display('tmp(5,1) = '+string(df51)+';');
display('tmp(5,2) = '+string(df52)+';');
display('tmp(5,3) = '+string(df53)+';');
display('tmp(5,4) = '+string(df54)+';');
display('tmp(5,5) = '+string(df55)+';');

%%
for i = 1: 5
    for j = 1:5
        display('df'+string(i)+string(j)+' = diff(df'+string(i)+', C'+string(j)+')');
    end
end