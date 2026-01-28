#!/usr/bin/env python3
"""
Biomni API Configuration Switcher

This script provides an easy way to switch between different API configuration profiles.
"""

import os
import sys
from pathlib import Path


class ProfileManager:
    """Manages API configuration profiles for Biomni."""

    def __init__(self, project_root: Path = None):
        """Initialize the profile manager.

        Args:
            project_root: Path to the Biomni project root. If None, uses current directory.
        """
        if project_root is None:
            # Assume we're in the Biomni directory
            self.project_root = Path.cwd()
        else:
            self.project_root = Path(project_root)

        self.profiles_dir = self.project_root / "profiles"
        self.env_file = self.project_root / ".env"

    def list_profiles(self) -> list[str]:
        """List all available profile files.

        Returns:
            List of profile names (without .env extension)
        """
        if not self.profiles_dir.exists():
            return []

        profiles = []
        for file in self.profiles_dir.glob("*.env"):
            # Skip README and default template
            if file.name not in ["README.md"]:
                profiles.append(file.stem)
        return sorted(profiles)

    def get_current_profile(self) -> str | None:
        """Get the name of the currently active profile.

        Returns:
            Profile name if currently using a profile, None otherwise
        """
        if not self.env_file.exists():
            return None

        # Check if .env is a symlink to a profile
        if self.env_file.is_symlink():
            target = self.env_file.resolve()
            if target.parent == self.profiles_dir:
                return target.stem

        # Try to identify by checking content
        if self.env_file.exists():
            content = self.env_file.read_text()
            for profile in self.list_profiles():
                profile_file = self.profiles_dir / f"{profile}.env"
                if profile_file.exists():
                    profile_content = profile_file.read_text()
                    # Simple comparison (ignoring whitespace differences)
                    if content.strip() == profile_content.strip():
                        return profile

        return None

    def switch_profile(self, profile_name: str, backup: bool = True) -> bool:
        """Switch to a specific profile.

        Args:
            profile_name: Name of the profile to switch to (without .env extension)
            backup: Whether to backup existing .env file before switching

        Returns:
            True if successful, False otherwise
        """
        profile_file = self.profiles_dir / f"{profile_name}.env"

        if not profile_file.exists():
            print(f"Error: Profile '{profile_name}' not found!")
            return False

        # Backup existing .env if requested
        if backup and self.env_file.exists():
            backup_file = self.project_root / ".env.backup"
            self.env_file.replace(backup_file)
            print(f"Backed up existing .env to .env.backup")

        # Copy profile to .env
        import shutil
        shutil.copy(profile_file, self.env_file)
        print(f"Switched to profile: {profile_name}")
        return True

    def show_profile_info(self, profile_name: str):
        """Show information about a profile.

        Args:
            profile_name: Name of the profile (without .env extension)
        """
        profile_file = self.profiles_dir / f"{profile_name}.env"

        if not profile_file.exists():
            print(f"Error: Profile '{profile_name}' not found!")
            return

        print(f"\n=== Profile: {profile_name} ===")
        print(f"File: {profile_file}")
        print("\nConfiguration preview:")
        print("-" * 60)

        content = profile_file.read_text()

        # Show non-comment lines that have content
        for line in content.split("\n"):
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                # Mask API keys for security
                key, value = line.split("=", 1)
                if "KEY" in key.upper() or "TOKEN" in key.upper():
                    if value and value != "":
                        value = "***" + value[-4:] if len(value) > 4 else "***"
                print(f"  {key}={value}")

        print("-" * 60)

    def interactive_menu(self):
        """Run an interactive menu for profile management."""
        while True:
            print("\n" + "=" * 60)
            print("Biomni API Configuration Manager")
            print("=" * 60)

            current = self.get_current_profile()
            print(f"\nCurrent profile: {current if current else 'None (using default or custom .env)'}")

            profiles = self.list_profiles()
            if not profiles:
                print("\nNo profiles found in profiles/ directory")
                return

            print("\nAvailable profiles:")
            for i, profile in enumerate(profiles, 1):
                marker = " [CURRENT]" if profile == current else ""
                print(f"  {i}. {profile}{marker}")

            print("\nOptions:")
            print("  <number> - Switch to profile")
            print("  i<number> - Show profile info (e.g., i1)")
            print("  c - Create new profile")
            print("  e - Edit .env file directly")
            print("  q - Quit")

            choice = input("\nYour choice: ").strip().lower()

            if choice == "q":
                print("Goodbye!")
                break
            elif choice == "c":
                self.create_profile_interactive()
            elif choice == "e":
                self.edit_env()
            elif choice.startswith("i") and choice[1:].isdigit():
                idx = int(choice[1:]) - 1
                if 0 <= idx < len(profiles):
                    self.show_profile_info(profiles[idx])
            elif choice.isdigit():
                idx = int(choice) - 1
                if 0 <= idx < len(profiles):
                    profile = profiles[idx]
                    if self.switch_profile(profile):
                        print(f"\nSuccess! Now using profile: {profile}")
                        print("\nYou can now run Biomni with this configuration.")
                else:
                    print("Invalid selection!")
            else:
                print("Invalid choice!")

    def create_profile_interactive(self):
        """Interactively create a new profile."""
        print("\n=== Create New Profile ===")
        name = input("Profile name (without .env extension): ").strip()

        if not name:
            print("Cancelled.")
            return

        profile_file = self.profiles_dir / f"{name}.env"
        if profile_file.exists():
            overwrite = input(f"Profile '{name}' already exists. Overwrite? (y/N): ").strip().lower()
            if overwrite != "y":
                print("Cancelled.")
                return

        # Ask which template to use
        templates = [p for p in self.list_profiles() if p not in ["default"]]
        print("\nChoose a template (or leave blank for empty profile):")
        for i, tpl in enumerate(templates, 1):
            print(f"  {i}. {tpl}")

        template_choice = input("Template number (or press Enter for empty): ").strip()

        if template_choice.isdigit() and 0 < int(template_choice) <= len(templates):
            template_file = self.profiles_dir / f"{templates[int(template_choice) - 1]}.env"
            content = template_file.read_text()
        else:
            content = "# Biomni Configuration Profile\n"

        # Create/edit the profile
        print(f"\nCreating profile: {name}")
        print("Enter your configuration (press Ctrl+D or Ctrl+Z when done):")

        try:
            lines = []
            while True:
                try:
                    line = input()
                    lines.append(line)
                except EOFError:
                    break

            if lines:
                content = "\n".join(lines)

            profile_file.write_text(content)
            print(f"\nProfile '{name}' created successfully!")

        except KeyboardInterrupt:
            print("\nCancelled.")

    def edit_env(self):
        """Edit the .env file directly."""
        import subprocess

        editor = os.environ.get("EDITOR", "nano")

        print(f"\nOpening .env with {editor}...")
        try:
            subprocess.call([editor, str(self.env_file)])
        except FileNotFoundError:
            print(f"Error: Editor '{editor}' not found. Please edit .env manually.")


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Biomni API Configuration Manager",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                    # Interactive mode
  %(prog)s list              # List all profiles
  %(prog)s switch anthropic  # Switch to anthropic profile
  %(prog)s info openai       # Show profile info
        """
    )

    parser.add_argument(
        "command",
        nargs="?",
        default="interactive",
        help="Command to run (list, switch, info, or leave blank for interactive)",
    )

    parser.add_argument(
        "profile",
        nargs="?",
        help="Profile name (for switch/info commands)",
    )

    args = parser.parse_args()

    manager = ProfileManager()

    if args.command == "list":
        profiles = manager.list_profiles()
        if not profiles:
            print("No profiles found.")
            return

        current = manager.get_current_profile()
        print("Available profiles:")
        for profile in profiles:
            marker = " [CURRENT]" if profile == current else ""
            print(f"  - {profile}{marker}")

    elif args.command == "switch":
        if not args.profile:
            print("Error: Please specify a profile name")
            print("Usage: python switch_profile.py switch <profile_name>")
            return

        if manager.switch_profile(args.profile):
            print(f"Success! Now using profile: {args.profile}")
        else:
            sys.exit(1)

    elif args.command == "info":
        if not args.profile:
            print("Error: Please specify a profile name")
            print("Usage: python switch_profile.py info <profile_name>")
            return

        manager.show_profile_info(args.profile)

    elif args.command == "interactive" or not args.command:
        manager.interactive_menu()

    else:
        print(f"Unknown command: {args.command}")
        print("Run: python switch_profile.py --help")
        sys.exit(1)


if __name__ == "__main__":
    main()
