function symBlockDecomposition_TEST_SUITE
%%

load('symBlockDecomposition_TEST_01', 'M0');
[rows, cols] = size(M0);
for i = 2 : cols
    if ( abs(M0(1, i)) > 0 ) 
        break;  
    end
end

symBlockDecomposition_TEST(M0, i - 1);

%%

load('symBlockDecomposition_TEST_02', 'M0');
[rows, cols] = size(M0);
for i = 2 : cols
    if ( abs(M0(1, i)) > 0 ) 
        break;  
    end
end

symBlockDecomposition_TEST(M0, i - 1);

%%

load('symBlockDecomposition_TEST_03', 'M0');
[rows, cols] = size(M0);
for i = 2 : cols
    if ( abs(M0(1, i)) > 0 ) 
        break;  
    end
end

symBlockDecomposition_TEST(M0, i - 1);

%%

load('symBlockDecomposition_TEST_04', 'M0');
[rows, cols] = size(M0);
for i = 2 : cols
    if ( abs(M0(1, i)) > 0 ) 
        break;  
    end
end

symBlockDecomposition_TEST(M0, i - 1);





%%
end

