function target_id = find_nearest_visible_agent(ax, ay, bx_list, by_list, range)
    % Find nearest visible agent from a list
    dists = sqrt((bx_list - ax).^2 + (by_list - ay).^2);
    visible_ids = find(dists <= range);
    if isempty(visible_ids)
        target_id = [];
    else
        [~, min_idx] = min(dists(visible_ids));
        target_id = visible_ids(min_idx);
    end
end
