function available = find_idle_ergodic(erg_flags)
    available = find(~erg_flags, 1); % pick first non-ergodic agent
end