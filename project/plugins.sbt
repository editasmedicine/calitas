resolvers += Resolver.url("fix-sbt-plugin-releases", url("https://dl.bintray.com/sbt/sbt-plugin-releases"))(Resolver.ivyStylePatterns)
addDependencyTreePlugin
addSbtPlugin("com.typesafe.sbt"  % "sbt-git"       % "1.0.0")
addSbtPlugin("com.eed3si9n"      % "sbt-assembly"  % "1.2.0")
addSbtPlugin("org.scoverage"     % "sbt-scoverage" % "1.6.1")