object FormGravity: TFormGravity
  Left = 322
  Top = 122
  AutoSize = True
  Caption = 'Gravity'
  ClientHeight = 169
  ClientWidth = 161
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object GBGravity: TGroupBox
    Left = 0
    Top = 0
    Width = 161
    Height = 169
    Caption = 'gravity vector'
    Color = clMoneyGreen
    ParentColor = False
    TabOrder = 0
    object Lgx: TLabel
      Left = 16
      Top = 32
      Width = 11
      Height = 13
      Caption = 'gx'
    end
    object Lgy: TLabel
      Left = 16
      Top = 64
      Width = 11
      Height = 13
      Caption = 'gy'
    end
    object Lgz: TLabel
      Left = 16
      Top = 96
      Width = 11
      Height = 13
      Caption = 'gz'
    end
    object lblgx: TLabel
      Left = 112
      Top = 24
      Width = 24
      Height = 13
      Caption = 'm/s2'
    end
    object lblgy: TLabel
      Left = 112
      Top = 56
      Width = 24
      Height = 13
      Caption = 'm/s2'
    end
    object lbldz: TLabel
      Left = 112
      Top = 88
      Width = 24
      Height = 13
      Caption = 'm/s2'
    end
    object Egx: TEdit
      Left = 40
      Top = 24
      Width = 65
      Height = 21
      TabOrder = 0
    end
    object Egy: TEdit
      Left = 40
      Top = 56
      Width = 65
      Height = 21
      TabOrder = 1
    end
    object Egz: TEdit
      Left = 40
      Top = 88
      Width = 65
      Height = 21
      TabOrder = 2
    end
    object Bapply: TButton
      Left = 56
      Top = 128
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 3
      OnClick = BapplyClick
    end
  end
end
