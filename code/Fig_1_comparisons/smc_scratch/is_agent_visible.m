function visible = is_agent_visible(ax, ay, bx, by, range)
    visible = norm([ax - bx, ay - by]) <= range;
end