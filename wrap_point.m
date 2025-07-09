function yn = wrap_point(blk, x, y)

    [~, r] = in_domain(blk, x, y);
    yn = y - r*blk.pitch;

end