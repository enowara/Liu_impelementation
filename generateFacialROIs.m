% generate ~ 15-sh landmark point regions
function [facialRegions] = generateFacialROIs(firstPoints)

p2_3 = [(firstPoints(2,1) + firstPoints(3,1))/2; (firstPoints(2,2) +firstPoints(3,2))/2];
%     half way between 2 and 3 
p32_40 = [(firstPoints(32,1) + firstPoints(40,1))/2; (firstPoints(32,2) +firstPoints(40,2))/2];
p36_43 = [(firstPoints(36,1) + firstPoints(43,1))/2; (firstPoints(36,2) +firstPoints(43,2))/2];

p_4_32 = [(firstPoints(4,1) + firstPoints(32,1))/2; (firstPoints(4,2) +firstPoints(32,2))/2];

p1_49a = [(p2_3(1)+p32_40(1))/2; (p2_3(2)+p32_40(2))/2];

p_15_16 = [(firstPoints(15,1) + firstPoints(16,1))/2; (firstPoints(15,2) +firstPoints(16,2))/2];

p17_35a = [(p36_43(1) + p_15_16(1))/2; (p36_43(2)+p_15_16(2))/2];
p17_35b = [(firstPoints(36,1) + firstPoints(14,1))/2; (firstPoints(36,2) +firstPoints(14,2))/2];



x1 = [firstPoints(1,1); p2_3(1); p1_49a(1)];
y1 = [firstPoints(1,2); p2_3(2); p1_49a(2)];
    region1b = boundary(x1,y1);
    region1 = [x1(region1b), y1(region1b)];
%region 2
x2 = [firstPoints(1,1); firstPoints(40,1); p32_40(1); p1_49a(1)];  
y2 = [firstPoints(1,2); firstPoints(40,2); p32_40(2); p1_49a(2)];  
    region2b = boundary(x2,y2);
    region2 = [x2(region2b), y2(region2b)];
% region 3
x3 = [firstPoints(4,1); p2_3(1); p_4_32(1); p1_49a(1)];
y3 = [firstPoints(4,2); p2_3(2); p_4_32(2); p1_49a(2)];
    region3b = boundary(x3,y3);
    region3 = [x3(region3b), y3(region3b)];
% region 4
x4 = [p1_49a(1); p_4_32(1); p32_40(1); firstPoints(32,1)];
y4 = [p1_49a(2); p_4_32(2); p32_40(2); firstPoints(32,2)];
    region4b = boundary(x4,y4);
    region4 = [x4(region4b), y4(region4b)];
% region 5

x5 = [p_4_32(1); firstPoints(4,1); firstPoints(6,1); firstPoints(49,1)];
y5 = [p_4_32(2); firstPoints(4,2); firstPoints(6,2); firstPoints(49,2)];
    region5b = boundary(x5,y5);
    region5 = [x5(region5b), y5(region5b)];
    
% region 6        
x6 = [p_4_32(1); firstPoints(32,1); firstPoints(50,1); firstPoints(49,1)];
y6 = [p_4_32(2); firstPoints(32,2); firstPoints(50,2); firstPoints(49,2)];
    region6b = boundary(x6,y6);
    region6 = [x6(region6b), y6(region6b)]; 
    
% region 7        
x7 = [firstPoints(6,1); firstPoints(49,1); firstPoints(58,1); firstPoints(9,1)];
y7 = [firstPoints(6,2); firstPoints(49,2); firstPoints(58,2); firstPoints(9,2)];
    region7b = boundary(x7,y7);
    region7 = [x7(region7b), y7(region7b)];     
    
% region 8        
x8 = [firstPoints(9,1); firstPoints(12,1); firstPoints(55,1); firstPoints(58,1)];
y8 = [firstPoints(9,2); firstPoints(12,2); firstPoints(55,2); firstPoints(58,2)];
    region8b = boundary(x8,y8);
    region8 = [x8(region8b), y8(region8b)];    
    
% region 9  FIX ME 
x9 = [firstPoints(12,1); firstPoints(14,1); firstPoints(55,1); p17_35b(1)];
y9 = [firstPoints(12,2); firstPoints(14,2); firstPoints(55,2); p17_35b(2)];
    region9b = boundary(x9,y9);
    region9 = [x9(region9b), y9(region9b)];    
    
% region 10
x10 = [firstPoints(14,1); p_15_16(1) ; p17_35a(1); p17_35b(1)];
y10 = [firstPoints(14,2); p_15_16(2) ; p17_35a(2); p17_35b(2)];
    region10b = boundary(x10,y10);
    region10 = [x10(region10b), y10(region10b)];     
    
% region 11
x11 = [firstPoints(55,1); firstPoints(54,1); firstPoints(36,1); p17_35b(1)];
y11 = [firstPoints(55,2); firstPoints(54,2); firstPoints(36,2); p17_35b(2)];
    region11b = boundary(x11,y11);
    region11 = [x11(region11b), y11(region11b)];   
    
% region 12    FIX ME
x12 = [p36_43(1); firstPoints(36,1); p17_35a(1); p17_35b(1)];
y12 = [p36_43(2); firstPoints(36,2); p17_35a(2); p17_35b(2)];
    region12b = boundary(x12,y12);
    region12 = [x11(region12b), y11(region12b)];  
    
% region 13    FIX ME 
x13 = [firstPoints(17,1); p17_35a(1); p_15_16(1)];
y13 = [firstPoints(17,2); p17_35a(2); p_15_16(2)];
    region13b = boundary(x13,y13);
    region13 = [x13(region13b), y13(region13b)];      
    
% region 14    FIX ME 
x14 = [firstPoints(43,1); firstPoints(17,1); p17_35a(1); p36_43(1) ];
y14 = [firstPoints(43,2); firstPoints(17,2); p17_35a(2); p36_43(2) ];
    region14b = boundary(x14,y14);
    region14 = [x14(region14b), y14(region14b)];      
    
    
% region 15    
x15 = [firstPoints(32,1);  firstPoints(36,1); firstPoints(40,1); firstPoints(43,1)];
y15 = [firstPoints(32,2);  firstPoints(36,2); firstPoints(40,2); firstPoints(43,2)];
    region15b = boundary(x15,y15);
    region15 = [x15(region15b), y15(region15b)];    
    
% region 16    
x16 = [firstPoints(32,1);  firstPoints(36,1);  firstPoints(55,1);  firstPoints(54,1); ...
    firstPoints(49,1);  firstPoints(50,1)]; 
y16 =  [firstPoints(32,2);  firstPoints(36,2);  firstPoints(55,2);  firstPoints(54,2); ...
    firstPoints(49,2);  firstPoints(50,2)]; 

    region16b = boundary(x16,y16);
    region16 = [x16(region16b), y16(region16b)];    


facialRegions{1} = region1;
facialRegions{2} = region2;
facialRegions{3} = region3;
facialRegions{4} = region4;
facialRegions{5} = region5;
facialRegions{6} = region6;
facialRegions{7} = region7;
facialRegions{8} = region8;
facialRegions{9} = region9;
facialRegions{10} = region10;
facialRegions{11} = region11;
facialRegions{12} = region12;
facialRegions{13} = region13;
facialRegions{14} = region14;
facialRegions{15} = region15;
facialRegions{16} = region16;


end
   
    
    