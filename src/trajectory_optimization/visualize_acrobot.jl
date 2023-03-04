function visualize_acrobot_trajectory(robot, xtraj, utraj, Thist)
    vis = Visualizer()
    open(vis)
    TrajOptPlots.set_mesh!(vis, robot)
    X1 = [SVector{4}(x) for x in eachcol(xtraj)];
    visualize!(vis, robot, Thist[end], X1)
end