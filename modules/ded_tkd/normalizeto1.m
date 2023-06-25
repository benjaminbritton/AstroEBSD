function normalized_matrix = normalizeto1(matrix_in)
normalized_matrix = matrix_in - min(matrix_in(:)); % make the lowest 0
normalized_matrix = normalized_matrix ./ max(normalized_matrix(:)); % Normalize
end
