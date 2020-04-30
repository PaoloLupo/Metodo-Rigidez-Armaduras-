
%% FASE DE ENTRADA DE DATOS %%
%datos de los nodos (x,y)
nodos = [
    0 0
    120 0
    120 120
    240 0
];

%elementos (i, j, E, A)
elementos = [
    1 3 30*10^6 3 
    1 2 30*10^6 3 
    2 3 30*10^6 3
    2 4 30*10^6 3
    4 3 30*10^6 3
];

% (nodo,(x=1, y=2),fuerza)
fuerzas = [
    2 1 0
    2 2 0
    3 1 8*10^3
    3 2 -8*10^3
    4 1 0
];



%% FASE SOLUCION %%


% matriz de fuerzas F
F = [];
F(:,1) = fuerzas(:,3)

% matrices ENSAMBLADA Y REDUCIDA
k_emsamblada = K_em(nodos, elementos)
k_reducida = K_red(nodos, elementos, fuerzas)


%matriz inversa de la reducida

k_invred = inv(k_reducida)

% multiplicacion de la inversa por la matriz de fuerzas
%u_des= desplazamientos desconocidos
F;

u_des = k_invred * F

nreduc = n_red(fuerzas);

%matriz columna de dezplazamientos
u = zeros(2*size(nodos,1),1);

for index = 1:size(fuerzas,1)
    u(nreduc(index)) = u_des(index);
    
end

u

%multiplicacion matriz ensamblada por los desplazamientos
%matriz de fuerzas final
F_des = k_emsamblada * u



%matriz de los nodos desplazados
nodos_des = nodos + transpose(reshape(u,2,size(nodos,1)))




%% GRAFICAS %%
% armadura inicial
x1 = nodos(:,1);
y1 = nodos(:,2);
x2 = [];
y2 = [];
for index = 1:size(elementos,1)

    x2 = [nodos(elementos(index,1),1) nodos(elementos(index,2),1)]; 
    y2 = [nodos(elementos(index,1),2) nodos(elementos(index,2),2)];
    p =plot(x2,y2,'-b'); hold on
    
end
q = plot(x1,y1,'*b');

%armadura modificada
x3 = nodos_des(:,1);
y3 = nodos_des(:,2);
x4 = [];
y4 = [];
for index = 1:size(elementos,1)

    x4 = [nodos_des(elementos(index,1),1) nodos_des(elementos(index,2),1)]; 
    y4 = [nodos_des(elementos(index,1),2) nodos_des(elementos(index,2),2)];
    r = plot(x4,y4,'-r'); hold on
    
end
s = plot(x3,y3,'*r');





%% FUNCIONES %%
% FASE MATRIZ DE RIGIDEZ EMSAMBLADA
%matriz de nodos por elemento
function xy = XY(nodos,elementos,N_elem)
    xy = [
        nodos(elementos(N_elem,1),:) 
        nodos(elementos(N_elem,2),:)
    ];
    
end

%matriz elementos 
function kel = KEL(XY, E, A)
    x1 = XY(1,1);
    y1 = XY(1,2);
    x2 = XY(2,1);
    y2 = XY(2,2);
    L = sqrt((x2-x1)^2 + (y2-y1)^2);
    kel = E*A/L * [1 0 -1 0; 0 0 0 0;-1 0 1 0; 0 0 0 0];

end

%matrices globales de cada elemento
function kgel = KGEL(XY,E,A)
    x1 = XY(1,1);
    y1 = XY(1,2);
    x2 = XY(2,1);
    y2 = XY(2,2);
    L = sqrt((x2-x1)^2 + (y2-y1)^2);
    if x2-x1 == 0
        c = 0;
        s = 1;
    else
        theta = atan((y2-y1)/(x2-x1));
        c = cos(theta);
        s = sin(theta);
    end
    kgel = E*A/L * [
        c^2 s*c -c^2 -s*c
        s*c s^2 -s*c -s^2
        -c^2 -s*c c^2 s*c
        -s*c -s^2 s*c s^2
                        
    ];
        
end

%funcion posicion para la matriz emsamblada
function gdl = GDL(elementos,N_elem)
    nodo_gdl_i = elementos(N_elem,1);
    nodo_gdl_j = elementos(N_elem,2);
    gdl = [
        2*nodo_gdl_i-1
        2*nodo_gdl_i
        2*nodo_gdl_j-1
        2*nodo_gdl_j
    ];
    
end

%funcion para hallar matriz a sumar en la emsamblada de cada elemento
function k = K_i(nodos, elementos, N_elem)
    k = zeros(2*length(nodos),2*length(nodos));
    k_i = KGEL(XY(nodos,elementos,N_elem), elementos(N_elem,3),elementos(N_elem,4));
    gdl_i = GDL(elementos,N_elem);
    
    for a = 1:4
        for b = 1:4
            k(gdl_i(a),gdl_i(b)) = k_i(a,b);
            
        end  
    end
end

% funcion matriz emsamblada de rigidez 
function k_em = K_em(nodos,elementos)
    k_em = zeros(2*length(nodos),2*length(nodos));
    for index = 1:size(elementos,1)
        ki = K_i(nodos, elementos, index);
        k_em = k_em + ki;
        
    end
end

%% FASE SOLUCIÒN

% funcion para hallar las filas y columnas para la reduccion de la matriz
function nred = n_red(fuerzas)
    nred = [];
    for index = 1:size(fuerzas,1)
        m = [
            fuerzas(index,1)*2-2 + fuerzas(index,2)
        ];
        nred(index) = m;
    end
end

% funcion para determinar la matriz reducida
function kred = K_red(nodos,elementos,fuerzas)
    k_emsamblada = K_em(nodos,elementos);
    for a = 1:size(fuerzas,1)
        for b = 1:size(fuerzas,1)
            z = n_red(fuerzas);
            w = z(a);
            y = z(b);
            kred(a,b) = k_emsamblada(w,y);
        end    
    end
    
end

