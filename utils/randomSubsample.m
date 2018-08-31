%--------------------------------------------------------------------------
function population = randomSubsample(population, n)
    % randomly select (with replacement) N objects from POPULATION
    if isempty(population)
        return;
    end

    if numel(population) > n
        population = population(ceil(numel(population) .* rand(n, 1))); % TW
    end
end % function
