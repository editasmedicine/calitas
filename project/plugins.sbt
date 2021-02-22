resolvers += Resolver.url("fix-sbt-plugin-releases", url("http://dl.bintray.com/sbt/sbt-plugin-releases"))(Resolver.ivyStylePatterns)

addSbtPlugin("com.typesafe.sbt"  % "sbt-git"       % "0.9.3")
addSbtPlugin("com.eed3si9n"      % "sbt-assembly"  % "0.14.6")
addSbtPlugin("org.scoverage"     % "sbt-scoverage" % "1.5.1")
addSbtPlugin("net.virtual-void"  % "sbt-dependency-graph" % "0.9.0")
