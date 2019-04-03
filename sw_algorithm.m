function [opt_score,a] = sw_algorithm(s1,s2,gp,scoring_m)
    bn_letters=['A';'R';'N';'D';'C';'Q';'E';'G';'H';'I';'L';'K';'N';'F';'P';'S';'T';'W';'Y';'V'];
    bn=scoring_m;
    
    % B funciona como um espaço no vetor dado que nao é uma letra presente
    % em BLOSUM
    s1 = ['B';s1];
    s2 = ['B';s2];
  
    l_1= length(s1);
    l_2= length(s2);
  
    % matriz com os scores
    sw_matrix = zeros(l_1,l_2);
  
    % matriz que indica as direções das setas
    % 1: seta para cima
    % 2: seta para a esquerda
    % 3: seta na diagonal
    d_matrix = zeros(l_1,l_2);
  
    d_vector=zeros(1,4);
  
    for i=2:l_1
        for j=2:l_2
            
            % procura indice da primeira letra na matriz BLOSUM
            for iter = 1:length(bn_letters)
                if(bn_letters(iter)==s1(i))
                    pos_s1=iter;
                    break;
                end
            end
            % procura indice da segunda letra na matriz BLOSUM
            for iter = 1:length(bn_letters)
                if(bn_letters(iter)==s2(j))
                    pos_s2=iter;
                    break;
                end
            end
            
            match=bn(pos_s1,pos_s2);
            
            % coloca em d_vector as diferentes possibilidades de score
            d_vector(1)=sw_matrix(i-1,j)+gp;
            d_vector(2)=sw_matrix(i,j-1)+gp;
            d_vector(3)=sw_matrix(i-1,j-1)+match;
            
            % coloca o maximo do vetor anterior como score definitivo na
            % matriz sw_matrix
            [M,I]=max(d_vector);
            sw_matrix(i,j)=M;
            
            % coloca em d_matrix um valor entre 1 e 3 correspondente à
            % direçao da seta, ou seja, de onde veio o score maximo
            if(M~=0)
                d_matrix(i,j)=I;
            end
        end
    end
    
    % encontra o maximo da matrix que corresponde ao score otimo
    % guarda em x e em y as coordenadas com score otimo
    opt_score = max(max(sw_matrix));
    [x,y]=find(sw_matrix==opt_score);
    x=x';
    y=y';
    j=1;
    % ciclo que cria tantos alinhamentos quantos os valores de score otimo
    for k=1:length(x)
        num=d_matrix(x(k),y(k));
        i=1;
    
        while num>0
            if(num==3)
                a(j,i)=s1(x(k));
                a(j+1,i)=s2(y(k));
                x(k)=x(k)-1;
                y(k)=y(k)-1;
            end
            if(num==1)
                a(j,i)=s1(x(k));
                a(j+1,i)='-';
                x(k)=x(k)-1;
            end
            if(num==2)
                a(j,i)='-';
                a(j+1,i)=s2(y(k));
                y(k)=y(k)-1;
            end
            num=d_matrix(x(k),y(k));
            i=i+1;
        end
        j=j+2;
    end
    
    a = fliplr(a);
end