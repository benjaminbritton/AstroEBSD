function normalized_matrix = normalizeto16bit(matrix_in)
% normalized_matrix = matrix_in;
normalized_matrix = matrix_in - min(matrix_in(:)); % make the lowest 0
normalized_matrix = normalized_matrix ./ max(normalized_matrix(:)); % Normalize
normalized_matrix = uint16(normalized_matrix * (2^16 - 1)); % Convert to 16 bit
end

