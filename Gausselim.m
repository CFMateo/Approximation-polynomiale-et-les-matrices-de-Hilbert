function x = Gausselim(A, b, i)
    U = A;  
    m = size(A, 1);  
    x = zeros(m, 1);
    
    for k = 1:m-1
        if i ~= 0  % Si i est différent de 0,  pivotage complet
            % On determine l'indice de la valeur la plus extreme dans la 
            % sous-matrice U(k:m, k:m):
            [max_val, max_index] = max(abs(U(k:m, k:m)), [], 'all', 'linear');
            [ind_row, ind_col] = ind2sub([m - k + 1, m - k + 1], max_index);
            
            % On ajuste les indices:
            ind_row = ind_row + k - 1;
            ind_col = ind_col + k - 1;
            
            % Échange de lignes pour amener le maximum à la position (k, k)
            U([k, ind_row], :) = U([ind_row, k], :);
            % On applique les mêmes échanges de lignes à b pour conserver
            % la cohérence:
            b([k, ind_row]) = b([ind_row, k]); 
            
            % Échange de colonnes pour amener le maximum
            % sur la diagonale principale (pivot)
            U(:, [k, ind_col]) = U(:, [ind_col, k]);
        end
        
        % Élimination de Gauss pour annuler les éléments en dessous du 
        % pivot dans la colonne k et transformer U en une matrice 
        % triangulaire supérieure:
        b(k+1:m) = b(k+1:m) - U(k+1:m, k) * b(k) / U(k, k);
        U(k+1:m, :) = U(k+1:m, :) - U(k+1:m, k) * U(k, :) / U(k, k);
    end
    
    % Résolution par substitution arrière pour trouver le vecteur solution 
    % x en partant de la dernière ligne jusqu'à la première:
    for j = m:-1:1
        x(j) = (b(j) - U(j, j+1:m) * x(j+1:m)) / U(j, j);
    end
end
