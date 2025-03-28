# --- Libraries ---
library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- Core Functions ---
custom_sign <- function(x) ifelse(x >= 0, 1, -1)

hopfield_update <- function(s, W) custom_sign(W %*% s)

converge_state <- function(s, W, max_iter = 100) {
  for (i in seq_len(max_iter)) {
    new_s <- hopfield_update(s, W)
    if (all(new_s == s)) break
    s <- new_s
  }
  s
}

enumerate_states <- function(N) as.matrix(expand.grid(rep(list(c(-1,1)), N)))

compute_basin_size <- function(W, target, N) {
  states <- enumerate_states(N)
  sum(apply(states, 1, function(s) all(converge_state(s,W)==target)))
}

get_basin_membership <- function(W, target, N) {
  states <- enumerate_states(N)
  labels <- apply(states, 1, paste, collapse="")
  labels[apply(states,1,function(s)all(converge_state(s,W)==target))]
}

# --- Simulation ---
run_single_simulation <- function(T2, N, p, p_pres) {
  memory_pres <- p_pres * matrix(sample(c(-1,1), p*N, replace=TRUE), nrow=p)
  repeat {
    target <- sample(c(-1,1), N, replace=TRUE)
    if (!any(apply(memory_pres,1,function(row) all(row==target)))) break
  }
  basin_sizes <- numeric(T2)
  W_list <- vector("list", T2)
  for (j in 1:T2) {
    memories <- rbind(memory_pres, matrix(rep(target,j), nrow=j, byrow=TRUE))
    W <- t(memories) %*% memories; diag(W) <- 0
    W_list[[j]] <- W
    basin_sizes[j] <- compute_basin_size(W,target,N)
  }
  list(basin_sizes=basin_sizes, W_list=W_list, target=target)
}

# --- Run & Detect Jumps ---
set.seed(1600)
N <- 8; T2 <- 1000; p <- 50; p_pres <- 10; threshold <- 10

sim_result <- run_single_simulation(T2,N,p,p_pres)
basin_sizes <- sim_result$basin_sizes
W_list <- sim_result$W_list
target <- sim_result$target

jump_events <- which(diff(basin_sizes) > threshold)
cat("Jump events detected at steps:",jump_events,"\n")

# User selects jump:
jump_chosen <- jump_events[1]
W_pre <- W_list[[jump_chosen]]
W_post <- W_list[[jump_chosen+1]]

# --- Graph Data ---
basin_pre <- get_basin_membership(W_pre,target,N)
basin_post <- get_basin_membership(W_post,target,N)

graph_from_W <- function(W, attractor_label, N){
  states <- enumerate_states(N)
  labels <- apply(states,1,paste,collapse="")
  edges <- data.frame(from=labels,
                      to=apply(states,1,function(s)paste(hopfield_update(s,W),collapse="")))
  g <- graph_from_data_frame(edges)
  g_sub <- induced_subgraph(g,subcomponent(g, attractor_label, mode="in"))
  reverse_edges(g_sub,E(g_sub))
}

label_target <- paste(target,collapse="")

g_pre <- graph_from_W(W_pre,label_target,N)
g_post <- graph_from_W(W_post,label_target,N)

# --- Unified layout (all nodes from both graphs) ---
unified_graph <- igraph::union(g_pre, g_post)

unified_layout <- create_layout(unified_graph, layout = "tree", circular = FALSE) %>%
  select(name, x, y) %>%
  mutate(y = max(y)-y) # root at bottom

vertices_pre <- data.frame(name=V(g_pre)$name) %>%
  left_join(unified_layout, by="name") %>%
  mutate(time="Pre", status="Old")

vertices_post <- data.frame(name=V(g_post)$name) %>%
  left_join(unified_layout, by="name") %>%
  mutate(time="Post", status=ifelse(name %in% basin_pre, "Old", "New"))

vertices_combined <- bind_rows(vertices_pre, vertices_post)

edges_pre <- igraph::as_data_frame(g_pre,"edges") %>% mutate(time="Pre",circular=FALSE,edge.id=row_number())
edges_post <- igraph::as_data_frame(g_post,"edges") %>% mutate(time="Post",circular=FALSE,edge.id=row_number())
edges_combined <- bind_rows(edges_pre,edges_post)

edges_combined <- edges_combined %>%
  left_join(vertices_combined %>% select(name,x,y,time),by=c("from"="name","time")) %>% rename(x=x,y=y) %>%
  left_join(vertices_combined %>% select(name,x,y,time),by=c("to"="name","time"),suffix=c("","end")) %>% rename(xend=xend,yend=yend)

# --- Plot ---
ggplot(vertices_combined,aes(x=x,y=y)) +
  geom_edge_diagonal(data=edges_combined,
                     aes(x=x,y=y,xend=xend,yend=yend,circular=circular,edge.id=edge.id),color="gray") +
  geom_point(aes(color=status),size=4) +
  facet_wrap(~factor(time,c("Pre","Post"))) +
  scale_color_manual(values=c(Old="forestgreen",New="dodgerblue")) +
  labs(title=paste("Target Attractor Basin from Step",jump_chosen,"to",jump_chosen+1),color="Status") +
  theme_minimal(base_size=14) +
  theme(plot.title=element_text(hjust=0.5,face="bold"))

ggplot(vertices_combined,aes(x=x,y=-y)) +
  geom_edge_diagonal(data=edges_combined,aes(x=x,y=-y,xend=xend,yend=-yend,circular=circular,edge.id=edge.id),color="gray")+
  geom_point(aes(color=status),size=4)+
  facet_wrap(~factor(time,c("Pre","Post")))+
  scale_color_manual(values=c(Old="forestgreen",New="dodgerblue"))+
  #labs(title=paste("Target Attractor Basin from Step",jump_chosen,"to",jump_chosen+1),color="Status")+
  theme_minimal(base_size=14)+
  theme(plot.title=element_text(hjust=0.5,face="bold"))

ggplot(vertices_combined, aes(x = x, y = y)) +
  geom_edge_diagonal(
    data = edges_combined,
    aes(x = x, y = y, xend = xend, yend = yend, circular = circular, edge.id = edge.id),
    color = "gray"
  ) +
  geom_point(aes(color = status), size = 4) +
  facet_wrap(~factor(time, c("Pre", "Post"))) +
  scale_color_manual(values = c(Old = "forestgreen", New = "dodgerblue")) +
  theme_minimal(base_size = 14) +
  labs(color="States")+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )



# Update facet labels
vertices_combined$time <- factor(vertices_combined$time,
                                 levels = c("Pre", "Post"),
                                 labels = c("Before Increase", "After Increase"))

edges_combined$time <- factor(edges_combined$time,
                              levels = c("Pre", "Post"),
                              labels = c("Before Increase", "After Increase"))

# Plot
plot3=ggplot(vertices_combined, aes(x = x, y = y)) +
  geom_edge_diagonal(
    data = edges_combined,
    aes(x = x, y = y, xend = xend, yend = yend, circular = circular, edge.id = edge.id),
    color = "gray"
  ) +
  geom_point(aes(color = status), size = 4) +
  facet_wrap(~time) +
  scale_color_manual(values = c(Old = "forestgreen", New = "dodgerblue")) +
  theme_minimal(base_size = 14) +
  labs(color="States")+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.background = element_rect(fill = 'white',color='white'),
    plot.background = element_rect(fill = 'white',color='white')
  )

ggsave("/Users/Jonathon/Documents/Papers/Hopfield/Fig3.png",units='cm',width=20,height=20,dpi=900,plot=plot3)


